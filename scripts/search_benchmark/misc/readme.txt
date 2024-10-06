mkdir -p nogit/large_data/id20/fas/

wget ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt -O nogit/large_data/id20/fas/cath-domain-list.txt

wget ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/non-redundant-data-sets/cath-dataset-nonredundant-S20.fa -O nogit/large_data/id20/fas/cath-dataset-nonredundant-S20.fa

# https://scop.berkeley.edu/astral/subsets/ver=2.08&seqOption=1
# から
# nogit/large_data/id20/fas/astral-scopedom-seqres-sel-gs-bib-20-2.08.fa
# を取得し保存。これはマニュアルで行った。

perl scripts/search_benchmark/misc/assign_cath.pl 
perl scripts/search_benchmark/misc/to_upper.pl 
/home/ubuntu8/apps/ncbi-blast-2.16.0+/bin/makeblastdb -dbtype prot -parse_seqids -in nogit/large_data/id20/fas/scope20_upper.fas
mkdir -p nogit/large_data/id20/fas/blast
mv nogit/large_data/id20/fas/scope20_upper.fas.* nogit/large_data/id20/fas/blast/
/home/ubuntu8/apps/ncbi-blast-2.16.0+/bin/blastp -query nogit/large_data/id20/fas/cath20_upper.fas -db nogit/large_data/id20/fas/blast/scope20_upper.fas -outfmt 6 -out nogit/large_data/id20/fas/blast/cath20_to_scope20.blastres.dat -num_threads 12
perl scripts/search_benchmark/misc/blastidfilter.pl

# を順に実行。
