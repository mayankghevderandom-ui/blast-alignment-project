
import pandas as pd

columns = [
    "Query", "Subject", "Identity", "Length",
    "Mismatches", "Gap openings", "Query start",
    "Query end", "Subject start", "Subject end",
    "E-value", "Bit score"
]

results = pd.read_csv(
    "results/blast_panel.txt",
    sep="\t",
    names=columns
)

print(results.to_string(index=False))
