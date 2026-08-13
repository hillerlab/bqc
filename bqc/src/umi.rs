// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! UMI extraction and relocation.
//!
//! Mirrors fastp's `--umi` feature: the UMI is read from a read prefix or a
//! stored index, attached to both mate names, and — for read-derived UMIs — the
//! prefix (plus any skip bases) is physically removed from the sequence.
//!
//! This is *not* UMI family clustering, error correction or consensus calling.
//! Those need biological grouping (alignment coordinates) and belong elsewhere.

use serde::Serialize;

use binseq::BinseqRecord;

use crate::error::{Error, Result};
use crate::read::ReadView;

/// Where a record's UMI lives.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
#[value(rename_all = "snake_case")]
pub enum UmiLocation {
    Index1,
    Index2,
    Read1,
    Read2,
    PerIndex,
    PerRead,
}

impl UmiLocation {
    /// Whether this location needs a second mate.
    #[must_use]
    pub fn needs_paired(self) -> bool {
        matches!(
            self,
            Self::Index2 | Self::Read2 | Self::PerIndex | Self::PerRead
        )
    }

    /// Whether this location physically removes sequence bases.
    #[must_use]
    pub fn removes_sequence(self) -> bool {
        matches!(self, Self::Read1 | Self::Read2 | Self::PerRead)
    }

    /// Stable name for reports and error messages.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Index1 => "index1",
            Self::Index2 => "index2",
            Self::Read1 => "read1",
            Self::Read2 => "read2",
            Self::PerIndex => "per_index",
            Self::PerRead => "per_read",
        }
    }
}
#[derive(Debug, Clone, Serialize)]
pub struct UmiStage {
    pub location: UmiLocation,
    pub length: usize,
    pub skip: usize,
    pub prefix: Vec<u8>,
    pub delimiter: Vec<u8>,
}

/// What UMI extraction did to one record's sequences.
#[derive(Debug, Clone, Copy)]
pub struct UmiOutcome {
    /// Bases removed from the front of R1 (0 for index modes).
    pub r1_clip: usize,
    /// Bases removed from the front of R2 (0 for SE/index modes).
    pub r2_clip: usize,
}

/// Worker-local UMI buffers, reused across records.
#[derive(Debug, Default)]
pub struct UmiScratch {
    /// The built tag: `delimiter` + optional `prefix_` + UMI.
    tag: Vec<u8>,
    /// R1 header with the tag inserted.
    r1_header: Vec<u8>,
    /// R2 header with the tag inserted.
    r2_header: Vec<u8>,
}

impl UmiScratch {
    /// The rewritten R1 header.
    #[must_use]
    pub fn r1_override(&self) -> &[u8] {
        &self.r1_header
    }

    /// The rewritten R2 header.
    #[must_use]
    pub fn r2_override(&self) -> &[u8] {
        &self.r2_header
    }
}

impl UmiStage {
    /// Extracts the UMI, computing per-mate clips and rewriting both headers.
    ///
    /// The rewritten headers live in `scratch` and stay valid until the next
    /// record is extracted. Sequence views are not mutated here: the caller
    /// clips them using the returned [`UmiOutcome`].
    pub fn extract<R: BinseqRecord>(
        &self,
        record: &R,
        r1: ReadView<'_>,
        r2: Option<ReadView<'_>>,
        scratch: &mut UmiScratch,
    ) -> Result<UmiOutcome> {
        let index = record.index();
        let (r1_clip, r2_clip) = match self.location {
            UmiLocation::Index1 => {
                let umi = first_index(record.sheader()).ok_or_else(|| {
                    Error::config(format!(
                        "record {index} has no parseable index in its R1 header"
                    ))
                })?;
                self.build_tag(&[umi], &mut scratch.tag);
                (0, 0)
            }
            UmiLocation::Index2 => {
                let umi = last_index(record.xheader()).ok_or_else(|| {
                    Error::config(format!(
                        "record {index} has no parseable index in its R2 header"
                    ))
                })?;
                self.build_tag(&[umi], &mut scratch.tag);
                (0, 0)
            }
            UmiLocation::Read1 => {
                let (clip, umi) = read_umi(r1.sequence, self.length, self.skip, index, "R1")?;
                self.build_tag(&[umi], &mut scratch.tag);
                (clip, 0)
            }
            UmiLocation::Read2 => {
                let r2 =
                    r2.ok_or_else(|| Error::config("read2 UMI requires a paired input file"))?;
                let (clip, umi) = read_umi(r2.sequence, self.length, self.skip, index, "R2")?;
                self.build_tag(&[umi], &mut scratch.tag);
                (0, clip)
            }
            UmiLocation::PerIndex => {
                let i1 = first_index(record.sheader()).ok_or_else(|| {
                    Error::config(format!(
                        "record {index} has no parseable index in its R1 header"
                    ))
                })?;
                let i2 = last_index(record.xheader()).ok_or_else(|| {
                    Error::config(format!(
                        "record {index} has no parseable index in its R2 header"
                    ))
                })?;
                self.build_tag(&[i1, i2], &mut scratch.tag);
                (0, 0)
            }
            UmiLocation::PerRead => {
                let r2 =
                    r2.ok_or_else(|| Error::config("per_read UMI requires a paired input file"))?;
                let (c1, u1) = read_umi(r1.sequence, self.length, self.skip, index, "R1")?;
                let (c2, u2) = read_umi(r2.sequence, self.length, self.skip, index, "R2")?;
                self.build_tag(&[u1, u2], &mut scratch.tag);
                (c1, c2)
            }
        };

        rewrite_header(record.sheader(), &scratch.tag, &mut scratch.r1_header);
        if r2.is_some() {
            rewrite_header(record.xheader(), &scratch.tag, &mut scratch.r2_header);
        }
        Ok(UmiOutcome { r1_clip, r2_clip })
    }

    /// Writes `delimiter` + optional `prefix_` + `parts` joined by `_` into `out`.
    fn build_tag(&self, parts: &[&[u8]], out: &mut Vec<u8>) {
        out.clear();
        out.extend_from_slice(&self.delimiter);
        if !self.prefix.is_empty() {
            out.extend_from_slice(&self.prefix);
            out.push(b'_');
        }
        for (i, part) in parts.iter().enumerate() {
            if i > 0 {
                out.push(b'_');
            }
            out.extend_from_slice(part);
        }
    }
}

/// Extracts the first `length` bases of `sequence` as the UMI, clipping
/// `length + skip` bases in total. Unlike fastp (which silently truncates a
/// short read's UMI), a read shorter than `length + skip` is an error.
fn read_umi<'a>(
    sequence: &'a [u8],
    length: usize,
    skip: usize,
    record: u64,
    mate: &'static str,
) -> Result<(usize, &'a [u8])> {
    let clip = length + skip;
    if sequence.len() < clip {
        return Err(Error::config(format!(
            "record {record} ({mate}) is {} bases long, shorter than the {length} bp UMI \
             plus {skip} skip base(s)",
            sequence.len()
        )));
    }
    Ok((clip, &sequence[..length]))
}

/// The index1 field: between the last `+` and the `:` before it (fastp's
/// `firstIndex`). Colon-only headers yield the field after the last `:`, the
/// same quirk fastp has.
fn first_index(header: &[u8]) -> Option<&[u8]> {
    if header.len() < 5 {
        return None;
    }
    let mut end = header.len();
    for i in (0..header.len() - 2).rev() {
        if header[i] == b'+' {
            end = i;
        }
        if header[i] == b':' {
            return Some(&header[i + 1..end]);
        }
    }
    None
}

/// The index2 field: everything after the last `:` or `+` (fastp's `lastIndex`).
fn last_index(header: &[u8]) -> Option<&[u8]> {
    if header.len() < 5 {
        return None;
    }
    let start = header.len() - 3;
    let pos = (0..=start)
        .rev()
        .find(|&i| header[i] == b':' || header[i] == b'+')?;
    Some(&header[pos + 1..])
}

/// Inserts `tag` before the first space in `header`, appending it when there is
/// no space, exactly as fastp does.
fn rewrite_header(header: &[u8], tag: &[u8], out: &mut Vec<u8>) {
    out.clear();
    if let Some(space) = header.iter().position(|&b| b == b' ') {
        out.extend_from_slice(&header[..space]);
        out.extend_from_slice(tag);
        out.extend_from_slice(&header[space..]);
    } else {
        out.extend_from_slice(header);
        out.extend_from_slice(tag);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_parsing_matches_fastp() {
        let header = b"NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA";
        assert_eq!(first_index(header), Some(&b"TATAGCCT"[..]));
        assert_eq!(last_index(header), Some(&b"GGTCCCGA"[..]));
    }

    #[test]
    fn header_rewrite_inserts_before_first_space_or_appends() {
        let mut out = Vec::new();
        rewrite_header(b"read1 1:N:0", b":ACGT", &mut out);
        assert_eq!(out, b"read1:ACGT 1:N:0");
        rewrite_header(b"read1", b":ACGT", &mut out);
        assert_eq!(out, b"read1:ACGT");
    }
}
