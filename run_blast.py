
import subprocess
from pathlib import Path

# Input and output paths
query_file = Path("data/query.fasta")
database_file = Path("data/protein_panel.fasta")
database_prefix = Path("data/protein_panel_db")
output_file = Path("results/blast_panel.txt")

# Check that the required input files exist
if not query_file.exists():
    raise FileNotFoundError(f"Query not found: {query_file}")

if not database_file.exists():
    raise FileNotFoundError(f"Protein panel not found: {database_file}")

# Create the results directory if necessary
output_file.parent.mkdir(parents=True, exist_ok=True)

# Build the local protein database
subprocess.run(
    [
        "makeblastdb",
        "-in", str(database_file),
        "-dbtype", "prot",
        "-out", str(database_prefix),
        "-parse_seqids",
    ],
    check=True,
)

# Search the database with BLASTP
subprocess.run(
    [
        "blastp",
        "-query", str(query_file),
        "-db", str(database_prefix),
        "-out", str(output_file),
        "-outfmt", "6",
        "-evalue", "1e-5",
    ],
    check=True,
)

print("Local BLAST search completed.")
print(f"Results saved to: {output_file}")
print("\nBLAST results:")

with output_file.open() as results:
    print(results.read(), end="")

