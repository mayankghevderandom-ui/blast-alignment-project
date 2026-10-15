
from Bio import Entrez, SeqIO

# Enter your email address
Entrez.email = "mayankghevde.random@gmail.com"

# Search NCBI for green fluorescent protein from Aequorea victoria
handle = Entrez.esearch(
    db="protein",
    term='"green fluorescent protein"[Protein Name] AND Aequorea victoria[Organism]',
    retmax=5
)

record = Entrez.read(handle)
handle.close()

ids = record["IdList"]

if not ids:
    raise RuntimeError("No GFP records found.")

# Download the first matching protein record
handle = Entrez.efetch(
    db="protein",
    id=ids[0],
    rettype="fasta",
    retmode="text"
)

protein = SeqIO.read(handle, "fasta")
handle.close()

# Save the protein sequence
SeqIO.write(protein, "data/source_GFP.fasta", "fasta")

print("Protein description:", protein.description)
print("Sequence length:", len(protein.seq), "amino acids")
print("Saved to: data/source_GFP.fasta")
print("Sequence:")
print(protein.seq)

