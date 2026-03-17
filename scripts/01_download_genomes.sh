#!/usr/bin/env bash

set -e

mkdir -p data/raw

for TYPE in 1 2 3 4
do

echo "Downloading DENV-$TYPE genomes..."

esearch -db nucleotide \
-query "\"Dengue virus $TYPE\"[Organism] AND complete genome" \
| efetch -format gb \
> data/raw/denv${TYPE}.gb

done

echo "Download finished"
