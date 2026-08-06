# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm AS builder

WORKDIR /app

COPY bqc/Cargo.toml bqc/Cargo.lock ./
COPY bqc/src ./src

RUN cargo build --release --all-features --bin bqc --locked && \
    strip target/release/bqc

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/bqc /usr/local/bin/bqc

# Set up non-root user
RUN useradd -m -u 1000 cuser && \
    chmod +x /usr/local/bin/bqc

USER cuser
WORKDIR /data

RUN bqc --help
