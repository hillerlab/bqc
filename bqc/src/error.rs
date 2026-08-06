// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Error types for `bqc`.
//!
//! The library reports every failure through [`Error`]. Variants are
//! user-oriented: each message is intended to be printed directly by the CLI
//! without additional context.

use std::fmt;
use std::path::{Path, PathBuf};

/// Convenience alias used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// All failures produced by `bqc`.
#[derive(Debug)]
pub enum Error {
    /// The input file is not a CBQ file.
    InvalidInputFormat { path: PathBuf, reason: String },
    /// A quality-dependent operation was requested on a quality-free file.
    MissingQuality(&'static str),
    /// A quality byte is outside the Phred+33 range.
    InvalidQualityEncoding {
        record: u64,
        mate: &'static str,
        byte: u8,
    },
    /// Sequence and quality lengths disagree for one mate.
    SequenceQualityLengthMismatch {
        record: u64,
        mate: &'static str,
        sequence: usize,
        quality: usize,
    },
    /// An adapter sequence is unusable.
    InvalidAdapter(String),
    /// The requested configuration is inconsistent, empty, or out of range.
    InvalidConfiguration(String),
    /// The output path already exists and `--force` was not supplied.
    OutputExists(PathBuf),
    /// Input and output refer to the same file.
    InputOutputConflict(PathBuf),
    /// A file could not be read.
    ReadError {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A file could not be written.
    WriteError {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A worker thread panicked.
    WorkerPanic,
    /// A failure reported by `binseq` itself.
    Binseq(binseq::Error),
    /// The CBQ input is structurally corrupt.
    CorruptInput { path: PathBuf, reason: String },
}

impl Error {
    /// Wraps an I/O failure that occurred while reading `path`.
    pub fn read(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Self::ReadError {
            path: path.as_ref().to_path_buf(),
            source,
        }
    }

    /// Wraps an I/O failure that occurred while writing `path`.
    pub fn write(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Self::WriteError {
            path: path.as_ref().to_path_buf(),
            source,
        }
    }

    /// Builds an [`Error::InvalidConfiguration`] from a formatted message.
    pub fn config(message: impl Into<String>) -> Self {
        Self::InvalidConfiguration(message.into())
    }

    /// Builds an [`Error::CorruptInput`] for `path`.
    pub fn corrupt(path: impl AsRef<Path>, reason: impl Into<String>) -> Self {
        Self::CorruptInput {
            path: path.as_ref().to_path_buf(),
            reason: reason.into(),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInputFormat { path, reason } => write!(
                f,
                "{} is not a readable CBQ file ({reason}); \
                 convert BQ/VBQ/FASTQ input with `bqtools encode` first",
                path.display()
            ),
            Self::MissingQuality(op) => write!(
                f,
                "quality-based operation requested ({op}), \
                 but input CBQ has no quality column"
            ),
            Self::InvalidQualityEncoding { record, mate, byte } => write!(
                f,
                "invalid quality byte {byte} in record {record} ({mate}): \
                 expected Phred+33 encoding (bytes 33..=126)"
            ),
            Self::SequenceQualityLengthMismatch {
                record,
                mate,
                sequence,
                quality,
            } => write!(
                f,
                "record {record} ({mate}) has {sequence} bases but {quality} quality values"
            ),
            Self::InvalidAdapter(reason) => write!(f, "invalid adapter: {reason}"),
            Self::InvalidConfiguration(reason) => write!(f, "{reason}"),
            Self::OutputExists(path) => write!(
                f,
                "refusing to overwrite existing file {} (use --force)",
                path.display()
            ),
            Self::InputOutputConflict(path) => {
                write!(f, "input and output are the same file: {}", path.display())
            }
            Self::ReadError { path, source } => {
                write!(f, "failed to read {}: {source}", path.display())
            }
            Self::WriteError { path, source } => {
                write!(f, "failed to write {}: {source}", path.display())
            }
            Self::WorkerPanic => write!(f, "a worker thread panicked; aborting"),
            Self::Binseq(source) => write!(f, "binseq error: {source}"),
            Self::CorruptInput { path, reason } => {
                write!(f, "corrupt CBQ file {}: {reason}", path.display())
            }
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReadError { source, .. } | Self::WriteError { source, .. } => Some(source),
            Self::Binseq(source) => Some(source),
            _ => None,
        }
    }
}

impl From<binseq::Error> for Error {
    fn from(source: binseq::Error) -> Self {
        Self::Binseq(source)
    }
}
