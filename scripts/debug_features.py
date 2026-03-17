from Bio import SeqIO

for record in SeqIO.parse("data/raw/denv1.gb", "genbank"):
    for feature in record.features:
        if feature.type == "CDS":
            print(feature.qualifiers)
            exit()
