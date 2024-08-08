# GMalign

Pointer to the state which had used in CASP16.

## Compile
```
# normal
cargo build --release

# with avx2
RUSTFLAGS="-C target-feature=+avx2" cargo build --release

# with sse3
RUSTFLAGS="-C target-feature=+sse3" cargo build --release
```

## Test
```
mkdir nogit # Some files will be created under this folder.

cargo test --release

# with avx2
RUSTFLAGS="-C target-feature=+avx2" cargo test --release

# with sse3
RUSTFLAGS="-C target-feature=+sse3" cargo test --release
```
