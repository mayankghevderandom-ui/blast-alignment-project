import random
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


query_file = Path("data/query.fasta")
shuffled_file = Path("data/shuffled_query.fasta")
database_prefix = Path("data/protein_panel_db")
output_file = Path("results/blast_shuffled.txt")
temporary_output = Path("results/.blast_shuffled_run.txt")

if not query_file.exists():
    raise FileNotFoundError(f"Query not found: {query_file}")

if not Path(str(database_prefix) + ".pin").exists():
    raise FileNotFoundError("Protein panel database not found. Run python run_blast.py first.")

original_record = SeqIO.read(query_file, "fasta")
original_sequence = str(original_record.seq)

output_file.parent.mkdir(parents=True, exist_ok=True)

total_shuffles = 100
shuffles_with_hits = 0

# Start with an empty results file.
output_file.write_text("")

for seed in range(total_shuffles):
    shuffled_amino_acids = list(original_sequence)
    random.Random(seed).shuffle(shuffled_amino_acids)
    shuffled_sequence = "".join(shuffled_amino_acids)

    # Confirm length and amino-acid composition are preserved.
    assert len(shuffled_sequence) == len(original_sequence)
    assert sorted(shuffled_sequence) == sorted(original_sequence)

    shuffled_record = SeqRecord(
        Seq(shuffled_sequence),
        id=f"shuffled_query_seed_{seed}",
        description="",
    )

    SeqIO.write(shuffled_record, shuffled_file, "fasta")

    subprocess.run(
        [
            "blastp",
            "-query", str(shuffled_file),
            "-db", str(database_prefix),
            "-out", str(temporary_output),
            "-outfmt", "6",
            "-evalue", "1e-5",
        ],
        check=True,
    )

    hits = temporary_output.read_text().strip()

    if hits:
        shuffles_with_hits += 1
        with output_file.open("a") as results:
            for hit in hits.splitlines():
                results.write(hit + "\n")

    temporary_output.unlink(missing_ok=True)

print(f"Original query length: {len(original_sequence)}")
print(f"Total shuffled sequences tested: {total_shuffles}")
print("Amino-acid composition preserved in every shuffle: True")
print(f"Shuffles producing a BLAST hit: {shuffles_with_hits}")
print(f"Shuffles producing no hits: {total_shuffles - shuffles_with_hits}")
print("E-value threshold: 1e-5")
print(f"Saved BLAST hit results: {output_file}")

if shuffles_with_hits == 0:
    print(f"Negative control passed: 0 of {total_shuffles} shuffles produced a hit.")
else:
    print("Review the reported hits before interpreting the negative control.")
