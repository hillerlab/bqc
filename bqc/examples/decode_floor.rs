// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Measures the decode-only floor of a CBQ file.
//!
//! ```bash
//! cargo run --release --example decode_floor -- input.cbq [threads]
//! ```
//!
//! This reads every record through exactly the path `bqc` uses — memory map,
//! block walk, zstd decompression, record iteration — and does no processing and
//! no writing. It is the lower bound on any run over that file, and therefore
//! the number that decides how much CPU-side optimization can ever be worth.

use std::hint::black_box;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use binseq::cbq::ColumnarBlock;
use binseq::BinseqRecord;
use bqc::io::CbqInput;

fn main() -> bqc::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: decode_floor <input.cbq> [threads]");
    let threads: usize = args.next().map_or(1, |value| value.parse().unwrap_or(1));

    let input = CbqInput::open(Path::new(&path))?;
    let schema = input.schema();
    let started = std::time::Instant::now();
    let next = AtomicUsize::new(0);
    let bases = AtomicUsize::new(0);
    let records = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| {
                let mut block = ColumnarBlock::new(input.header());
                let mut dctx = zstd::zstd_safe::DCtx::create();
                let mut local_bases = 0usize;
                let mut local_records = 0usize;
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    if index >= input.blocks().len() {
                        break;
                    }
                    let range = input.load(index, &mut block, &mut dctx).expect("decode");
                    for record in block.iter_records(range) {
                        // Touch every column the processing stages would read.
                        local_bases += black_box(record.sseq()).len();
                        if schema.paired {
                            local_bases += black_box(record.xseq()).len();
                        }
                        if schema.quality {
                            black_box(record.squal());
                            if schema.paired {
                                black_box(record.xqual());
                            }
                        }
                        if schema.headers {
                            black_box(record.sheader());
                        }
                        local_records += 1;
                    }
                }
                bases.fetch_add(local_bases, Ordering::Relaxed);
                records.fetch_add(local_records, Ordering::Relaxed);
            });
        }
    });

    let elapsed = started.elapsed().as_secs_f64();
    let records = records.load(Ordering::Relaxed);
    let bases = bases.load(Ordering::Relaxed);
    println!(
        "decode-only  T={threads}  {elapsed:.3} s  {:.0} records/s  {:.0} bases/s  \
         ({records} records, {bases} bases)",
        records as f64 / elapsed,
        bases as f64 / elapsed,
    );
    Ok(())
}
