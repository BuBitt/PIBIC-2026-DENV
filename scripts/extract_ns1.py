from Bio import SeqIO
import sys

if len(sys.argv) != 3:
    print("Usage: python extract_ns1.py input.gb output.fasta")
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

count = 0

with open(output_file, "w") as out:

    for record in SeqIO.parse(input_file, "genbank"):

        for feature in record.features:

            if feature.type in ["CDS", "mat_peptide"]:

                product = feature.qualifiers.get("product", [""])[0].lower()

                if "ns1" in product:

                    seq = feature.extract(record.seq)

                    out.write(f">{record.id}\n{seq}\n")

                    count += 1

print("NS1 sequences extracted:", count)
