# GPSM-align (General Position Specific * Matrix aligner)
Create multiple sequence alignment (MSA) with some position specific values of proteins. e.g. sequence representaion of protein language models, position specific scoring matrix, hmm profiles, and so on.

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

## Example
You have to install ESM2 https://github.com/facebookresearch/esm before running these examples. & have to download esm2_t33_650M_UR50D.pt

### Prepare Representations
```
mkdir -p nogit/run_example

python scripts/esm_process.py --infile example_files/esm2_650m_example_output/example_seq.fas --outdir nogit/run_example/esmout --crop_length 600 --shift_length 200 --cut_length 100  --model_path ＜path_to_the_esm2_t33_650M_UR50D.pt＞ --device cuda

# create a list of representations except for the first sequence. We will use this as "query" and list as "template DB".
find nogit/run_example/esmout|sort|grep mat.gz|tail -n 9 > nogit/run_example/esmout_list.dat
```

### Generate A3M Alignment (For Template Search For Template Based Modeling)
```
target/release/matrix_align --in nogit/run_example/esmout/seq_0.mat.gz --in_list nogit/run_example/esmout_list.dat --out nogit/run_example/example_res.a3m  --num_threads 32 --a3m_pairwise true --alignment_type local --score_type dot_product --normalize true
```
A3m alignment file will be saved in nogit/run_example/example_res.a3m.
With this command, score/profile_length > 300 may be homologous proteins.

### Generate MSA
```
target/release/matrix_align --in nogit/run_example/esmout/seq_0.mat.gz,nogit/run_example/esmout/seq_1.mat.gz,nogit/run_example/esmout/seq_2.mat.gz,nogit/run_example/esmout/seq_3.mat.gz --out nogit/run_example/example_res.msa  --num_threads 32 --a3m_pairwise false --alignment_type global --score_type dot_product --normalize true --tree_type NJ --distance_base averaged_value
```
MSA will be saved in nogit/run_example/example_res.msa.
With this command, "score" does not represent something much.
