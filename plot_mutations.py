
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# Read the reference sequence and mutation table
reference = SeqIO.read("data/source_GFP.fasta", "fasta")
mutations = pd.read_csv("data/mutations.csv")

sequence_length = len(reference.seq)

# Create the figure
fig, ax = plt.subplots(figsize=(12, 3))

# Draw the protein sequence as a horizontal line
ax.hlines(
    y=0,
    xmin=1,
    xmax=sequence_length,
    linewidth=4
)

# Mark and label each mutation
for _, row in mutations.iterrows():
    position = int(row["Position"])
    original = row["Original"]
    mutated = row["Mutated"]

    label = f"{original}{position}{mutated}"

    ax.scatter(position, 0, s=100, zorder=3)
    ax.annotate(
        label,
        xy=(position, 0),
        xytext=(0, 18),
        textcoords="offset points",
        ha="center"
    )

# Format the graph
ax.set_xlim(0, sequence_length + 1)
ax.set_ylim(-0.5, 0.8)
ax.set_yticks([])
ax.set_xlabel("Amino-acid position")
ax.set_title("Artificially Introduced Mutations in GFP")

plt.tight_layout()

# Save the figure
plt.savefig("results/mutation_positions.png", dpi=300)
plt.show()

print("Saved: results/mutation_positions.png")
