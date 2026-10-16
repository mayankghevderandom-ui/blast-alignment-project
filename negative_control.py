
import random
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Input and output files
query_file = Path("data/query.fasta")
shuffled_file = Path("data/shuffled_query.fasta")
database_prefix = Path("data/protein_panel_db")
output_file = Path("results/blast_shuffled.txt")

# Check that the required files exist
if not query_file.exists():
    raise FileNotFoundError(f"Query not found: {query_file}")

if not Path(str(database_prefix) + ".pin").exists():
    raise FileNotFoundError(
        "Protein panel database not found. Run python run_blast.py first."
    )

# Read the original query sequence
record = SeqIO.read(query_file, "fasta")
original_sequence = str(record.seq)

# Shuffle amino-acid order reproducibly
shuffled_amino_acids = list(original_sequence)
random.Random(42).shuffle(shuffled_amino_acids)
shuffled_sequence = "".join(shuffled_amino_acids)

# Save the shuffled sequence
shuffled_record = SeqRecord(
    Seq(shuffled_sequence),
    id="shuffled_query",
    description="Composition-preserving shuffled GFP query"
)
SeqIO.write(shuffled_record, shuffled_file, "fasta")

# Confirm length and amino-acid composition are preserved
assert len(shuffled_sequence) == len(original_sequence)
assert sorted(shuffled_sequence) == sorted(original_sequence)

# Run BLASTP against the same five-protein database
output_file.parent.mkdir(parents=True, exist_ok=True)

subprocess.run(
    [
        "blastp",
        "-query", str(shuffled_file),
        "-db", str(database_prefix),
        "-out", str(output_file),
        "-outfmt", "6",
        "-evalue", "1e-5",
    ],
    check=True,
)

# Report the results
print("Original query length:", len(original_sequence))
print("Shuffled query length:", len(shuffled_sequence))
print("Amino-acid composition preserved: True")
print("Saved shuffled sequence:", shuffled_file)
print("Saved BLAST results:", output_file)

if output_file.stat().st_size == 0:
    print("No hits reported at E-value threshold 1e-5.")
else:
    print("\nReported BLAST hits:")
    print(output_file.read_text(), end="")
