// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! `bqc` — CBQ-native adapter removal, read trimming and per-read filtering.
//!
//! The library is organised around one idea: a transformation is a pair of
//! coordinates into a borrowed record. [`read::Span`] carries those coordinates,
//! [`process::Workflow`] updates them, and [`io`] writes the corresponding
//! slices back out. Nothing in the hot path allocates or rewrites sequence and
//! quality data.
//!
//! # Guarantees
//!
//! * **Order preservation.** Output records appear in input order for any
//!   thread count. See [`engine`].
//! * **Schema preservation.** Pairing, qualities, headers and flags are taken
//!   from the input file header and reproduced exactly. Synthetic headers and
//!   quality values that `binseq` exposes for files without those columns are
//!   never written out.
//! * **Determinism.** Adapter tie-breaking, failure-reason ordering and stage
//!   order are fixed; the same input and configuration always produce the same
//!   decoded output.
//! * **Atomic outputs.** Every file is written to a temporary path and renamed
//!   only after a successful run.
//!
//! # Example
//!
//! ```no_run
//! use bqc::adapter::{Adapter, AdapterParams, AdapterStage};
//! use bqc::process::Workflow;
//! use bqc::read::ReadView;
//!
//! let stage = AdapterStage::new(
//!     vec![Adapter::new("illumina", b"AGATCGGAAGAGC")?],
//!     Vec::new(),
//!     AdapterParams::default(),
//!     None,
//! )?;
//! let workflow = Workflow::new(Some(stage), None, None, None, None)?;
//!
//! let sequence = b"ACGTACGTACGTACGTAGATCGGAAGAGC";
//! let read = ReadView::new(sequence, None, 0, "R1")?;
//! let result = workflow.process(0, read, None)?;
//! assert_eq!(result.r1.final_length(), 16);
//! # Ok::<(), bqc::error::Error>(())
//! ```

pub mod adapter;
pub mod cli;
pub mod config;
pub mod correct;
pub mod detect;
pub mod engine;
pub mod error;
pub mod filter;
pub mod io;
pub mod linked;
pub mod overlap;
pub mod process;
pub mod read;
pub mod report;
pub mod segment;
pub mod sniff;
pub mod stats;
pub mod trim;

pub use error::{Error, Result};
