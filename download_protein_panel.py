
from pathlib import Path
from urllib.request import Request, urlopen
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO

# Output file for all proteins
output_dir = Path("data")
output_dir.mkdir(exist_ok=True)

output_file = output_dir / "protein_panel.fasta"

# UniProt accession and descriptive name
proteins = {
    "P42212": "Aequorea_victoria_GFP",
    "Q9U6Y8": "Discosoma_red_fluorescent_protein",
    "P61626": "Human_lysozyme",
    "P01308": "Human_insulin",
    "P02144": "Human_myoglobin",
}

records = []

for accession, protein_name in proteins.items():
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = Request(
        url,
        headers={"User-Agent": "BLAST-Alignment-Mini-Project"}
    )

    print(f"Downloading {protein_name} ({accession})...")

    with urlopen(request, timeout=30) as response:
        fasta_text = response.read().decode("utf-8")

    record = SeqIO.read(StringIO(fasta_text), "fasta")
    record.id = accession
    record.description = protein_name
    records.append(record)

# Save all downloaded sequences in one multi-FASTA file
SeqIO.write(records, output_file, "fasta")

print(f"\nDownloaded {len(records)} proteins.")
print(f"Saved to: {output_file}")

print("\nSequence summary:")
for record in records:
    print(f"{record.id}\t{record.description}\t{len(record.seq)} amino acids")
