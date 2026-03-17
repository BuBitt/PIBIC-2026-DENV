#!/usr/bin/env bash

set -e

for TYPE in 1 2 3 4
do

echo "Processing DENV-$TYPE"

# baixar genomas
esearch -db nucleotide \
-query "\"Dengue virus $TYPE\"[Organism] AND complete genome" \
| efetch -format gb \
> data/raw/denv${TYPE}.gb


# extrair NS1
python scripts/extract_ns1.py \
data/raw/denv${TYPE}.gb \
data/intermediate/denv${TYPE}_NS1.fasta


# remover duplicatas
seqkit rmdup -s \
data/intermediate/denv${TYPE}_NS1.fasta \
> data/intermediate/denv${TYPE}_NS1_nodup.fasta


# cluster 99%
cd-hit-est \
-i data/intermediate/denv${TYPE}_NS1_nodup.fasta \
-o data/intermediate/denv${TYPE}_NS1_cluster99.fasta \
-c 0.99


# alinhamento
mafft --auto --thread -1 \
data/intermediate/denv${TYPE}_NS1_cluster99.fasta \
> data/aligned/denv${TYPE}_NS1_aligned.fasta


# filogenia
iqtree \
-s data/aligned/denv${TYPE}_NS1_aligned.fasta \
-m MFP \
-bb 1000 \
-nt AUTO \
-pre results/trees/denv${TYPE}

done
