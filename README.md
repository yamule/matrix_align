# GMalign

## Compile
```
# normal
cargo build --release

# with avx2
RUSTFLAGS="-C target-feature=+avx2" cargo build --release

# with sse3
RUSTFLAGS="-C target-feature=+sse3" cargo build --release
```

