from Bio import AlignIO
import numpy as np
import pandas as pd

# --- carregar alinhamento ---
alignment = AlignIO.read("denv4_NS1_aligned.fasta", "fasta")

# --- calcular conservação por posição ---
def conservation(column):
    counts = {}
    for aa in column:
        if aa != '-' and aa != 'X':
            counts[aa] = counts.get(aa, 0) + 1
    if not counts:
        return np.nan
    return max(counts.values()) / sum(counts.values())

scores = []
for i in range(alignment.get_alignment_length()):
    col = alignment[:, i]
    scores.append(conservation(col))

# --- epítopos (IEDB) ---
epitopes = [
    (26, 39),
    (62, 62),
    (93, 113),
    (116, 130),
    (138, 148),
    (172, 179),
    (210, 210),
    (228, 241),
    (248, 272),
    (291, 318),
    (338, 348),
]

# --- calcular conservação média dos epítopos ---
results = []
for start, end in epitopes:
    region = scores[start-1:end]
    mean_score = np.nanmean(region)
    results.append([start, end, mean_score])

df = pd.DataFrame(results, columns=["Start", "End", "Conservation"])

print(df)
df.to_csv("epitopos_conservacao_denv4.csv", index=False)
