
import random
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Make the random choices reproducible
random.seed(42)

# Load the original GFP sequence
source = SeqIO.read("data/source_GFP.fasta", "fasta")
original = str(source.seq)

# List of the 20 standard amino acid codes
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Select three distinct positions
positions = sorted(random.sample(range(len(original)), 3))

# Make a copy of the sequence as a list
query = list(original)
mutations = []

# Replace the amino acid at each selected position
for position in positions:
    old_amino_acid = original[position]

    possible_amino_acids = [
        aa for aa in amino_acids if aa != old_amino_acid
    ]

    new_amino_acid = random.choice(possible_amino_acids)
    query[position] = new_amino_acid

    mutations.append({
        "Position": position + 1,
        "Original": old_amino_acid,
        "Mutated": new_amino_acid
    })

# Convert the modified list back into a sequence string
query_sequence = "".join(query)

# Create a FASTA record for the unknown protein
unknown = SeqRecord(
    Seq(query_sequence),
    id="unknown_query",
    description="Artificially mutated GFP"
)

# Save the query sequence
SeqIO.write(unknown, "data/query.fasta", "fasta")

# Save the known mutations in a CSV table
mutation_table = pd.DataFrame(mutations)
mutation_table.to_csv("data/mutations.csv", index=False)

print("Original length:", len(original))
print("Query length:", len(query_sequence))
print("\nIntroduced mutations:")
print(mutation_table.to_string(index=False))
print("\nSaved: data/query.fasta")
print("Saved: data/mutations.csv")
