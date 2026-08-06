// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

use std::process::ExitCode;

use bqc::cli::{run, Cli, Outcome};
use clap::Parser;

/// Exit statuses, so a pipeline can branch without parsing stderr.
///
/// ```text
/// 0  the command completed
/// 2  a command line, configuration or runtime error
/// 3  the result was not confident, and --require-confident was given
/// ```
///
/// An inconclusive result is a successful analysis: without
/// `--require-confident` it exits 0.
const EXIT_ERROR: u8 = 2;
const EXIT_NOT_CONFIDENT: u8 = 3;

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(Outcome::Success) => ExitCode::SUCCESS,
        Ok(Outcome::NotConfident) => ExitCode::from(EXIT_NOT_CONFIDENT),
        Err(error) => {
            eprintln!("error: {error}");
            let mut source = std::error::Error::source(&error);
            while let Some(cause) = source {
                eprintln!("  caused by: {cause}");
                source = cause.source();
            }
            ExitCode::from(EXIT_ERROR)
        }
    }
}
