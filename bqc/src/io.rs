// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! CBQ input, schema-preserving output, atomic commits and the reason sidecar.
//!
//! Input blocks are located by walking the file's block headers, which makes
//! each CBQ block independently addressable. That is the unit of work used by
//! the processing engine: because block boundaries come from the input, the
//! output is partitioned identically regardless of thread count.

use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use binseq::cbq::{
    BlockHeader, BlockRange, ColumnarBlock, ColumnarBlockWriter, FileHeader, FileHeaderBuilder,
    IndexFooter, IndexHeader,
};
use binseq::{BinseqRecord, SequencingRecordBuilder};
use memmap2::Mmap;
use serde::Serialize;
use zstd::zstd_safe;

use crate::error::{Error, Result};
use crate::process::Mate;
use crate::read::Span;

const FILE_HEADER_SIZE: usize = size_of::<FileHeader>();
const BLOCK_HEADER_SIZE: usize = size_of::<BlockHeader>();
const INDEX_HEADER_SIZE: usize = size_of::<IndexHeader>();
const INDEX_FOOTER_SIZE: usize = size_of::<IndexFooter>();

/// An 8-byte aligned scratch buffer for `bytemuck`-backed header casts.
#[repr(align(8))]
struct Aligned<const N: usize>([u8; N]);

impl<const N: usize> Aligned<N> {
    fn from_slice(bytes: &[u8]) -> Self {
        let mut buf = Self([0u8; N]);
        buf.0.copy_from_slice(bytes);
        buf
    }
}

/// The data-presence attributes of a CBQ file.
///
/// The schema always comes from the file header, never from record contents:
/// `binseq` synthesizes headers and quality values for files that do not store
/// them, and those synthetic values must never be written back out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Schema {
    pub paired: bool,
    pub quality: bool,
    pub headers: bool,
    pub flags: bool,
}

impl Schema {
    #[must_use]
    pub fn from_header(header: &FileHeader) -> Self {
        Self {
            paired: header.is_paired(),
            quality: header.has_qualities(),
            headers: header.has_headers(),
            flags: header.has_flags(),
        }
    }

    /// Builds an output file header preserving this schema.
    #[must_use]
    pub fn to_file_header(self, compression_level: usize, block_size: usize) -> FileHeader {
        FileHeaderBuilder::default()
            .is_paired(self.paired)
            .with_qualities(self.quality)
            .with_headers(self.headers)
            .with_flags(self.flags)
            .with_compression_level(compression_level)
            .with_block_size(block_size)
            .build()
    }

    /// The single-end view of this schema, used for orphan outputs.
    #[must_use]
    pub fn unpaired(self) -> Self {
        Self {
            paired: false,
            ..self
        }
    }
}

/// One addressable CBQ block.
#[derive(Debug, Clone, Copy)]
pub struct InputBlock {
    header: BlockHeader,
    data_start: usize,
    data_end: usize,
    /// Global index of this block's first record.
    pub first_record: u64,
    pub num_records: u64,
}

impl InputBlock {
    /// Global index just past this block's last record.
    #[must_use]
    pub fn end_record(&self) -> u64 {
        self.first_record + self.num_records
    }
}

/// A memory-mapped CBQ input file.
pub struct CbqInput {
    path: PathBuf,
    mmap: Mmap,
    header: FileHeader,
    blocks: Vec<InputBlock>,
    num_records: u64,
}

impl CbqInput {
    /// Opens and indexes a CBQ file.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|e| Error::read(path, e))?;
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| Error::read(path, e))?;
        if mmap.len() < FILE_HEADER_SIZE {
            return Err(Error::InvalidInputFormat {
                path: path.to_path_buf(),
                reason: "file is shorter than a CBQ file header".to_string(),
            });
        }
        let header = FileHeader::from_bytes(
            &Aligned::<FILE_HEADER_SIZE>::from_slice(&mmap[..FILE_HEADER_SIZE]).0,
        )
        .map_err(|e| Error::InvalidInputFormat {
            path: path.to_path_buf(),
            reason: e.to_string(),
        })?;

        // Walk the block headers. The block stream ends at the index header,
        // whose magic differs, so a failed header parse marks the boundary.
        let mut blocks = Vec::new();
        let mut offset = FILE_HEADER_SIZE;
        let mut first_record = 0u64;
        while offset + BLOCK_HEADER_SIZE <= mmap.len() {
            let Ok(block_header) = BlockHeader::from_bytes(
                &Aligned::<BLOCK_HEADER_SIZE>::from_slice(
                    &mmap[offset..offset + BLOCK_HEADER_SIZE],
                )
                .0,
            ) else {
                break;
            };
            let data_start = offset + BLOCK_HEADER_SIZE;
            let data_end = data_start
                .checked_add(block_header.block_len())
                .filter(|end| *end <= mmap.len())
                .ok_or_else(|| {
                    Error::corrupt(
                        path,
                        format!("block at byte {offset} extends past end of file"),
                    )
                })?;
            blocks.push(InputBlock {
                header: block_header,
                data_start,
                data_end,
                first_record,
                num_records: block_header.num_records,
            });
            first_record = first_record
                .checked_add(block_header.num_records)
                .ok_or_else(|| {
                    Error::corrupt(path, "block record counts overflow the file record count")
                })?;
            offset = data_end;
        }

        // The block stream must be followed by a well-formed index trailer:
        // an `IndexHeader`, the compressed index and an `IndexFooter`. Without
        // this check a file truncated at a block boundary would be silently
        // accepted with only its leading records.
        let truncated = || Error::corrupt(path, "missing or truncated CBQ index");
        let index_header_end = offset
            .checked_add(INDEX_HEADER_SIZE)
            .ok_or_else(truncated)?;
        if index_header_end > mmap.len() {
            return Err(truncated());
        }
        IndexHeader::from_bytes(
            &Aligned::<INDEX_HEADER_SIZE>::from_slice(&mmap[offset..index_header_end]).0,
        )
        .map_err(|_| truncated())?;
        // `IndexHeader` keeps its sizes crate-private; read the compressed
        // length directly (CBQ format version 1 layout: magic, u_bytes, z_bytes).
        let z_bytes = u64::from_le_bytes(
            mmap[offset + 16..offset + 24]
                .try_into()
                .map_err(|_| truncated())?,
        ) as usize;
        let footer_start = index_header_end
            .checked_add(z_bytes)
            .ok_or_else(truncated)?;
        let file_end = footer_start
            .checked_add(INDEX_FOOTER_SIZE)
            .ok_or_else(truncated)?;
        if file_end != mmap.len() {
            return Err(truncated());
        }
        IndexFooter::from_bytes(
            &Aligned::<INDEX_FOOTER_SIZE>::from_slice(&mmap[footer_start..file_end]).0,
        )
        .map_err(|_| truncated())?;

        Ok(Self {
            path: path.to_path_buf(),
            mmap,
            header,
            blocks,
            num_records: first_record,
        })
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub fn header(&self) -> FileHeader {
        self.header
    }

    #[must_use]
    pub fn schema(&self) -> Schema {
        Schema::from_header(&self.header)
    }

    #[must_use]
    pub fn num_records(&self) -> u64 {
        self.num_records
    }

    #[must_use]
    pub fn blocks(&self) -> &[InputBlock] {
        &self.blocks
    }

    /// Decompresses one block into `block`, returning the range needed to
    /// iterate its records with correct global indices.
    pub fn load(
        &self,
        index: usize,
        block: &mut ColumnarBlock,
        dctx: &mut zstd_safe::DCtx<'_>,
    ) -> Result<BlockRange> {
        let input = self.blocks[index];
        block.decompress_from_bytes(
            &self.mmap[input.data_start..input.data_end],
            input.header,
            dctx,
        )?;
        Ok(BlockRange::new(0, input.end_record()))
    }
}

/// A worker-local buffer of finished, already-compressed output blocks.
pub type Fragment = ColumnarBlockWriter<Vec<u8>>;

/// Creates an output fragment that shares the destination file's schema.
pub fn fragment(header: FileHeader) -> Result<Fragment> {
    Ok(ColumnarBlockWriter::new_headless(Vec::new(), header)?)
}

/// What one mate contributes to an output record.
///
/// Normally just a span into the record's own bytes. Base correction replaces a
/// mate's sequence and quality with worker-owned buffers, and those cover the
/// whole mate, so the span indexes them exactly as it indexes the record.
#[derive(Debug, Clone, Copy)]
pub struct MateOutput<'a> {
    /// Span into the record's own bytes.
    pub span: Span,
    /// Corrected sequence, quality and the span into them. Present only when
    /// correction rewrote this mate; the span may use a different coordinate
    /// system than `span` after UMI clipping.
    corrected: Option<(&'a [u8], &'a [u8], Span)>,
    /// Header override, applied in place of the record's own header.
    pub header: Option<&'a [u8]>,
}

impl<'a> MateOutput<'a> {
    /// A mate written straight from the record.
    #[must_use]
    pub fn borrowed(span: Span) -> Self {
        Self {
            span,
            corrected: None,
            header: None,
        }
    }

    /// A mate whose bases were corrected, written from `sequence` and `quality`.
    fn corrected(record_span: Span, scratch_span: Span, sequence: &'a [u8], quality: &'a [u8]) -> Self {
        Self {
            span: record_span,
            corrected: Some((sequence, quality, scratch_span)),
            header: None,
        }
    }

    /// A mate written from correction scratch when it was corrected, and from the
    /// record otherwise. `record_span` indexes the record's bytes; `scratch_span`
    /// indexes the corrected buffers (the two differ only after UMI clipping).
    #[must_use]
    pub fn from_scratch(
        record_span: Span,
        scratch_span: Span,
        scratch: &'a crate::correct::CorrectionScratch,
        mate: Mate,
    ) -> Self {
        match scratch.corrected(mate) {
            Some((sequence, quality)) => {
                Self::corrected(record_span, scratch_span, sequence, quality)
            }
            None => Self::borrowed(record_span),
        }
    }

    /// Attaches a header override.
    #[must_use]
    pub fn with_header(mut self, header: Option<&'a [u8]>) -> Self {
        self.header = header;
        self
    }

    /// The sequence to write, the corrected quality when there is one, and the
    /// span indexing both.
    fn bytes(self, from_record: &'a [u8]) -> (&'a [u8], Option<&'a [u8]>, Span) {
        match self.corrected {
            Some((sequence, quality, span)) => (sequence, Some(quality), span),
            None => (from_record, None, self.span),
        }
    }
}

/// Rejects a sequence/quality pair whose lengths disagree.
fn check_len(index: u64, mate: &'static str, seq: &[u8], qual: &[u8]) -> Result<()> {
    if qual.len() != seq.len() {
        return Err(Error::SequenceQualityLengthMismatch {
            record: index,
            mate,
            sequence: seq.len(),
            quality: qual.len(),
        });
    }
    Ok(())
}

/// Writes one record to `fragment`, retaining only the given spans.
///
/// Headers, flags and pairing are propagated unchanged; quality values are
/// sliced with exactly the same coordinates as the sequence.
pub fn push_record<R: BinseqRecord>(
    fragment: &mut Fragment,
    schema: Schema,
    record: &R,
    r1: MateOutput<'_>,
    r2: MateOutput<'_>,
) -> Result<()> {
    let (s_seq, s_qual_override, s_span) = r1.bytes(record.sseq());
    let (x_seq, x_qual_override, x_span) = r2.bytes(record.xseq());
    let s_qual = if schema.quality {
        let qual = s_qual_override.unwrap_or_else(|| record.squal());
        check_len(record.index(), "R1", s_seq, qual)?;
        Some(&qual[s_span.start..s_span.end])
    } else {
        None
    };
    let x_qual = if schema.quality && schema.paired {
        let qual = x_qual_override.unwrap_or_else(|| record.xqual());
        check_len(record.index(), "R2", x_seq, qual)?;
        Some(&qual[x_span.start..x_span.end])
    } else {
        None
    };

    let builder = SequencingRecordBuilder::default()
        .s_seq(&s_seq[s_span.start..s_span.end])
        .opt_s_qual(s_qual)
        .opt_s_header(if schema.headers {
            Some(r1.header.unwrap_or_else(|| record.sheader()))
        } else {
            None
        })
        .opt_x_seq(if schema.paired {
            Some(&x_seq[x_span.start..x_span.end])
        } else {
            None
        })
        .opt_x_qual(x_qual)
        .opt_x_header(if schema.headers && schema.paired {
            Some(r2.header.unwrap_or_else(|| record.xheader()))
        } else {
            None
        })
        .opt_flag(if schema.flags { record.flag() } else { None });

    fragment.push(builder.build()?)?;
    Ok(())
}

/// Writes one mate of `record` to `fragment` as a single-end record.
///
/// This is the orphan-output path: the surviving mate of a broken pair keeps
/// its own sequence, quality, header and the record's flag, but the
/// destination schema is single-end.
pub fn push_mate<R: BinseqRecord>(
    fragment: &mut Fragment,
    schema: Schema,
    record: &R,
    mate: crate::process::Mate,
    output: MateOutput<'_>,
) -> Result<()> {
    debug_assert!(!schema.paired, "orphan outputs are single-end");
    let header_override = output.header;
    let (sequence, quality, header) = match mate {
        crate::process::Mate::R1 => (record.sseq(), record.squal(), record.sheader()),
        crate::process::Mate::R2 => (record.xseq(), record.xqual(), record.xheader()),
    };
    let (sequence, corrected_quality, span) = output.bytes(sequence);
    let quality = corrected_quality.unwrap_or(quality);
    let quality = if schema.quality {
        check_len(record.index(), mate.name(), sequence, quality)?;
        Some(&quality[span.start..span.end])
    } else {
        None
    };
    let builder = SequencingRecordBuilder::default()
        .s_seq(&sequence[span.start..span.end])
        .opt_s_qual(quality)
        .opt_s_header(if schema.headers {
            Some(header_override.unwrap_or(header))
        } else {
            None
        })
        .opt_flag(if schema.flags { record.flag() } else { None });
    fragment.push(builder.build()?)?;
    Ok(())
}

/// Writes one fragment of a segmented read as a single-end record.
///
/// `header` is the fragment's own header — the source header plus the provenance
/// suffix — and is written only when the schema carries headers; a header-free
/// input stays header-free, and the sidecar is the provenance in that case.
pub fn push_segment<R: BinseqRecord>(
    fragment: &mut Fragment,
    schema: Schema,
    record: &R,
    span: Span,
    header: &[u8],
) -> Result<()> {
    debug_assert!(!schema.paired, "segmentation output is single-end");
    let sequence = record.sseq();
    let quality = if schema.quality {
        let quality = record.squal();
        check_len(record.index(), "R1", sequence, quality)?;
        Some(&quality[span.start..span.end])
    } else {
        None
    };
    let builder = SequencingRecordBuilder::default()
        .s_seq(&sequence[span.start..span.end])
        .opt_s_qual(quality)
        .opt_s_header(if schema.headers { Some(header) } else { None })
        // Flags describe the source molecule, so every fragment inherits them.
        .opt_flag(if schema.flags { record.flag() } else { None });
    fragment.push(builder.build()?)?;
    Ok(())
}

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Atomically reserves an invocation-unique temporary file beside `path`.
fn create_temp(path: &Path) -> Result<(PathBuf, File)> {
    let name = path
        .file_name()
        .ok_or_else(|| Error::config(format!("{} is not a valid output path", path.display())))?;
    for _ in 0..128 {
        let sequence = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let tmp = path.with_file_name(format!(
            ".{}.bqc-tmp-{}-{sequence}",
            name.to_string_lossy(),
            std::process::id()
        ));
        match OpenOptions::new().write(true).create_new(true).open(&tmp) {
            Ok(file) => return Ok((tmp, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => return Err(Error::write(&tmp, error)),
        }
    }
    Err(Error::write(
        path,
        std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "could not reserve a unique temporary output file",
        ),
    ))
}

/// Rejects output paths that already exist or collide with the input.
fn check_output_path(path: &Path, input: &Path, force: bool) -> Result<()> {
    if same_file(path, input) {
        return Err(Error::InputOutputConflict(path.to_path_buf()));
    }
    if path.exists() && !force {
        return Err(Error::OutputExists(path.to_path_buf()));
    }
    Ok(())
}

/// Publishes a completed temporary file.
///
/// A hard link is the no-clobber commit point: because the temporary file sits
/// beside its destination, both names are on the same filesystem, and creating
/// the destination either succeeds atomically or reports that it already
/// exists. Once the link succeeds the output is durable under `path`; removing
/// the temporary name is cleanup and cannot invalidate that commit.
fn publish(tmp: &Path, path: &Path, force: bool) -> Result<()> {
    if force {
        return fs::rename(tmp, path).map_err(|error| Error::write(path, error));
    }
    match fs::hard_link(tmp, path) {
        Ok(()) => {
            let _ = fs::remove_file(tmp);
            Ok(())
        }
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            Err(Error::OutputExists(path.to_path_buf()))
        }
        Err(error) => Err(Error::write(path, error)),
    }
}

/// Whether two paths resolve to the same file.
fn same_file(a: &Path, b: &Path) -> bool {
    match (fs::canonicalize(a), fs::canonicalize(b)) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    }
}

/// A CBQ output file written atomically.
///
/// Blocks are appended in commit order; the temporary file is published only
/// after [`CbqOutput::finish`] succeeds. If the run fails, the temporary file is
/// removed on drop and no partial output is left behind.
pub struct CbqOutput {
    path: PathBuf,
    tmp: PathBuf,
    writer: ColumnarBlockWriter<File>,
    force: bool,
}

impl CbqOutput {
    /// Creates the output file, refusing to clobber existing data.
    pub fn create(path: &Path, input: &Path, header: FileHeader, force: bool) -> Result<Self> {
        check_output_path(path, input, force)?;
        let (tmp, file) = create_temp(path)?;
        let writer = ColumnarBlockWriter::new(file, header)?;
        Ok(Self {
            path: path.to_path_buf(),
            tmp,
            writer,
            force,
        })
    }

    /// Appends every completed block of `fragment` to the output.
    pub fn commit(&mut self, fragment: &mut Fragment) -> Result<()> {
        self.writer.ingest_completed(fragment)?;
        Ok(())
    }

    /// Writes the index and publishes the temporary file.
    pub fn finish(&mut self) -> Result<()> {
        self.writer.finish()?;
        publish(&self.tmp, &self.path, self.force)
    }

    /// The file header this output was created with. Worker fragments must use
    /// the same header so their blocks are byte-compatible with the file.
    #[must_use]
    pub fn header(&self) -> FileHeader {
        self.writer.header()
    }
}

impl Drop for CbqOutput {
    fn drop(&mut self) {
        // Also retries cleanup after a successful hard-link publication whose
        // first unlink attempt failed.
        let _ = fs::remove_file(&self.tmp);
    }
}

/// A text file written atomically, used for the reason sidecar and reports.
pub struct TextOutput {
    path: PathBuf,
    tmp: PathBuf,
    writer: BufWriter<File>,
    force: bool,
}

impl TextOutput {
    /// Creates the file, refusing to clobber existing data.
    pub fn create(path: &Path, input: &Path, force: bool) -> Result<Self> {
        check_output_path(path, input, force)?;
        let (tmp, file) = create_temp(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            tmp,
            writer: BufWriter::new(file),
            force,
        })
    }

    /// Appends raw bytes.
    pub fn write(&mut self, bytes: &[u8]) -> Result<()> {
        self.writer
            .write_all(bytes)
            .map_err(|e| Error::write(&self.path, e))
    }

    /// Flushes and publishes the temporary file.
    pub fn finish(&mut self) -> Result<()> {
        self.writer
            .flush()
            .map_err(|e| Error::write(&self.path, e))?;
        publish(&self.tmp, &self.path, self.force)
    }
}

impl Drop for TextOutput {
    fn drop(&mut self) {
        // Also retries cleanup after a successful hard-link publication whose
        // first unlink attempt failed.
        let _ = fs::remove_file(&self.tmp);
    }
}

/// Column header of the failure-reason sidecar.
pub const REASON_HEADER: &str = "record_index\tmate\tstatus\treasons\toriginal_length\t\
adapter_trimmed_length\tfinal_length\tadapter_name\tadapter_start\n";

/// Column header of the segmentation provenance sidecar.
pub const SEGMENT_HEADER: &str = "source_record_index\tsegment_index\tsource_mate\tstart\tend\t\
length\tleft_adapter\tright_adapter\toriginal_header\tstatus\tfilter_reasons\n";

/// Appends one provenance row, describing where a fragment came from.
///
/// Every emitted fragment gets a row, accepted or not: the sidecar is the only
/// record of provenance for a header-free input, so it cannot be partial.
pub fn write_segment_row<R: BinseqRecord>(
    out: &mut Vec<u8>,
    schema: Schema,
    record: &R,
    result: &crate::process::FragmentResult,
    adapter_name: impl Fn(usize) -> String,
) {
    use std::io::Write as _;
    let fragment = result.fragment;
    let name = |adapter: Option<usize>| match adapter {
        Some(index) => adapter_name(index),
        None => ".".to_string(),
    };
    let _ = write!(
        out,
        "{}\t{}\tR1\t{}\t{}\t{}\t{}\t{}\t",
        record.index(),
        fragment.index,
        result.retained.start,
        result.retained.end,
        result.retained.len(),
        name(fragment.left_adapter),
        name(fragment.right_adapter),
    );
    if schema.headers {
        write_escaped(out, record.sheader());
    } else {
        out.push(b'.');
    }
    let _ = writeln!(
        out,
        "\t{}\t{}",
        if result.passed() { "PASS" } else { "FAIL" },
        result.reasons.label(),
    );
}

/// Column header of the correction log, for the configured detail level.
#[must_use]
pub fn correction_log_header(detail: crate::correct::LogDetail) -> &'static str {
    match detail {
        crate::correct::LogDetail::Reads => {
            "record_index\tr1_header\tr2_header\toverlap_offset\toverlap_length\t\
             overlap_mismatches\tcorrected_r1_bases\tcorrected_r2_bases\t\
             unresolved_mismatches\tfinal_disposition\n"
        }
        crate::correct::LogDetail::Bases => {
            "record_index\tmate\tread_position\tdonor_mate\tdonor_position\t\
             original_base\tcorrected_base\toriginal_phred\tcorrected_phred\t\
             overlap_offset\toverlap_length\tfinal_disposition\n"
        }
    }
}

/// Writes `header` with tabs, newlines, carriage returns and backslashes escaped,
/// so one record can never spill into another row or column.
pub fn write_escaped(out: &mut Vec<u8>, header: &[u8]) {
    for &byte in header {
        match byte {
            b'\\' => out.extend_from_slice(b"\\\\"),
            b'\t' => out.extend_from_slice(b"\\t"),
            b'\n' => out.extend_from_slice(b"\\n"),
            b'\r' => out.extend_from_slice(b"\\r"),
            other => out.push(other),
        }
    }
}

/// Appends one sidecar row.
pub fn write_reason_row(
    out: &mut Vec<u8>,
    record_index: u64,
    mate: &str,
    result: &crate::process::MateResult,
    adapter_name: Option<&str>,
) {
    use std::io::Write as _;
    let status = if result.passed() { "PASS" } else { "FAIL" };
    let (name, start) = match (adapter_name, result.adapter_hit) {
        (Some(name), Some(hit)) => (name.to_string(), hit.start.to_string()),
        _ => (".".to_string(), ".".to_string()),
    };
    let _ = writeln!(
        out,
        "{record_index}\t{mate}\t{status}\t{}\t{}\t{}\t{}\t{name}\t{start}",
        result.reasons.label(),
        result.original_length,
        result.adapter_trimmed_length,
        result.final_length(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_sizes_match_the_cbq_layout() {
        assert_eq!(FILE_HEADER_SIZE, 64);
        assert_eq!(BLOCK_HEADER_SIZE, 96);
    }

    #[test]
    fn schema_round_trips_through_a_file_header() {
        for paired in [false, true] {
            for quality in [false, true] {
                for headers in [false, true] {
                    for flags in [false, true] {
                        let schema = Schema {
                            paired,
                            quality,
                            headers,
                            flags,
                        };
                        let header = schema.to_file_header(3, 4096);
                        assert_eq!(Schema::from_header(&header), schema);
                        assert_eq!(header.compression_level, 3);
                        assert_eq!(header.block_size, 4096);
                    }
                }
            }
        }
    }

    #[test]
    fn temporary_files_are_unique_and_stay_beside_the_target() {
        let directory = tempfile::tempdir().unwrap();
        let target = directory.path().join("clean.cbq");
        let (first, _first_file) = create_temp(&target).unwrap();
        let (second, _second_file) = create_temp(&target).unwrap();
        assert_eq!(first.parent(), Some(directory.path()));
        assert_eq!(second.parent(), Some(directory.path()));
        assert_ne!(first, second);
        assert!(
            first
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(".clean.cbq.bqc-tmp-")
        );
        assert!(create_temp(Path::new("/")).is_err());
    }

    #[test]
    fn no_clobber_publication_has_exactly_one_concurrent_winner() {
        use std::sync::{Arc, Barrier};

        let directory = tempfile::tempdir().unwrap();
        let target = directory.path().join("result.txt");
        let (first_tmp, mut first_file) = create_temp(&target).unwrap();
        let (second_tmp, mut second_file) = create_temp(&target).unwrap();
        first_file.write_all(b"first").unwrap();
        second_file.write_all(b"second").unwrap();
        drop(first_file);
        drop(second_file);

        let barrier = Arc::new(Barrier::new(3));
        let (first_result, second_result) = std::thread::scope(|scope| {
            let first_barrier = Arc::clone(&barrier);
            let first_target = target.clone();
            let first_tmp_for_thread = first_tmp.clone();
            let first = scope.spawn(move || {
                first_barrier.wait();
                publish(&first_tmp_for_thread, &first_target, false)
            });

            let second_barrier = Arc::clone(&barrier);
            let second_target = target.clone();
            let second_tmp_for_thread = second_tmp.clone();
            let second = scope.spawn(move || {
                second_barrier.wait();
                publish(&second_tmp_for_thread, &second_target, false)
            });

            barrier.wait();
            (first.join().unwrap(), second.join().unwrap())
        });

        let results = [&first_result, &second_result];
        assert_eq!(results.iter().filter(|result| result.is_ok()).count(), 1);
        assert_eq!(
            results
                .iter()
                .filter(
                    |result| matches!(result, Err(Error::OutputExists(path)) if path == &target)
                )
                .count(),
            1
        );
        let bytes = fs::read(&target).unwrap();
        assert!(bytes == b"first" || bytes == b"second");
        let _ = fs::remove_file(first_tmp);
        let _ = fs::remove_file(second_tmp);
    }

    #[test]
    fn forced_publication_replaces_the_destination() {
        let directory = tempfile::tempdir().unwrap();
        let target = directory.path().join("result.txt");
        fs::write(&target, b"old").unwrap();
        let (tmp, mut file) = create_temp(&target).unwrap();
        file.write_all(b"new").unwrap();
        drop(file);

        publish(&tmp, &target, true).unwrap();

        assert_eq!(fs::read(&target).unwrap(), b"new");
        assert!(!tmp.exists());
    }

    #[test]
    fn reason_rows_render_all_columns() {
        use crate::adapter::AdapterHit;
        use crate::filter::FilterReason;
        use crate::process::MateResult;
        use crate::trim::TrimOutcome;

        let result = MateResult {
            retained: Span { start: 0, end: 23 },
            umi_clip: 0,
            original_length: 151,
            adapter_trimmed_length: 40,
            adapter_hit: Some(AdapterHit {
                adapter: 0,
                start: 40,
                errors: 0,
                overlap: 33,
                consumed: 33,
            }),
            linked: None,
            overlap_removed: 0,
            trimming: TrimOutcome::default(),
            reasons: FilterReason::TOO_SHORT | FilterReason::TOO_MANY_N,
        };
        let mut out = Vec::new();
        write_reason_row(&mut out, 10291, "R1", &result, Some("illumina"));
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "10291\tR1\tFAIL\tTOO_SHORT/TOO_MANY_N\t151\t40\t23\tillumina\t40\n"
        );
    }

    #[test]
    fn passing_rows_report_pass_and_empty_adapter_columns() {
        use crate::filter::FilterReason;
        use crate::process::MateResult;
        use crate::trim::TrimOutcome;

        let result = MateResult {
            retained: Span { start: 0, end: 147 },
            umi_clip: 0,
            original_length: 151,
            adapter_trimmed_length: 151,
            adapter_hit: None,
            linked: None,
            overlap_removed: 0,
            trimming: TrimOutcome::default(),
            reasons: FilterReason::empty(),
        };
        let mut out = Vec::new();
        write_reason_row(&mut out, 10291, "R2", &result, None);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "10291\tR2\tPASS\tPASS\t151\t151\t147\t.\t.\n"
        );
    }
}
