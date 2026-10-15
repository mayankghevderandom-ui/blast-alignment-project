
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# Read the original and mutated protein sequences
reference = SeqIO.read("data/source_GFP.fasta", "fasta")
query = SeqIO.read("data/query.fasta", "fasta")

ref_seq = str(reference.seq)
query_seq = str(query.seq)

# Configure the pairwise aligner
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

# Align the sequences
alignment = aligner.align(ref_seq, query_seq)[0]

print("Alignment score:", alignment.score)
print(alignment)

# Identify substitutions
mismatches = []

for i, (ref_aa, query_aa) in enumerate(
    zip(ref_seq, query_seq), start=1
):
    if ref_aa != query_aa:
        mismatches.append({
            "Position": i,
            "Reference": ref_aa,
            "Query": query_aa,
            "Mutation": f"{ref_aa}{i}{query_aa}"
        })

# Save the mismatch table
table = pd.DataFrame(mismatches)
table.to_csv("results/mismatches.csv", index=False)

print("\nRecovered substitutions:")
print(table.to_string(index=False))

# Validate against the known mutation record
expected = pd.read_csv("data/mutations.csv")

found = set(zip(table["Position"], table["Query"]))
truth = set(zip(expected["Position"], expected["Mutated"]))

print("\nAll mutations recovered:", found == truth)
print("Saved: results/mismatches.csv")
