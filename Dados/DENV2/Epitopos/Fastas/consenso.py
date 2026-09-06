from Bio import AlignIO
import numpy as np
import pandas as pd

# --- carregar alinhamento ---
alignment = AlignIO.read("denv2_NS1_alinhado.fasta", "fasta")

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
    (10, 10),
    (26, 39),
    (93, 120),
    (122, 129),
    (140, 148),
    (172, 179),
    (227, 243),
    (248, 272),
    (278, 278),
    (280, 281),
    (290, 318),
    (324, 324),
    (338, 347),
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
