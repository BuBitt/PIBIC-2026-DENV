from Bio import AlignIO
import numpy as np
import pandas as pd

# --- carregar alinhamento ---
alignment = AlignIO.read("denv3_NS1_alinhado.fasta", "fasta")

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
    (26, 40),
    (80, 82),
    (93, 129),
    (141, 150),
    (173, 179),
    (228, 242),
    (249, 282),
    (290, 328),
    (337, 348),

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
