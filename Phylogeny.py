!pip install biopython

# Primate tree from cytochrome b protein sequences
# Question: Can protein sequences show how closely related animals are?

import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

Entrez.email = "memoonaanwar115@gmail.com"

# Animals we want to compare
animals = [
    "Homo sapiens",
    "Pan troglodytes",
    "Gorilla gorilla",
    "Pongo abelii",
    "Macaca mulatta",
    "Mus musculus",
    "Bos taurus",
]

# Step 1: download one protein for each animal
sequences = {}

for animal in animals:
    print("Getting", animal)
    search = animal + "[Organism] AND cytochrome b[Title] AND 370:390[SLEN]"

    handle = Entrez.esearch(db="protein", term=search, retmax=1)
    result = Entrez.read(handle)
    handle.close()

    if len(result["IdList"]) == 0:
        print("  not found, skipping")
        continue

    handle = Entrez.efetch(db="protein", id=result["IdList"][0],
                           rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()

    sequences[animal] = str(record.seq)
    print("  length:", len(record.seq))

# Step 2: cut each protein into pieces of 3 letters
def get_pieces(protein):
    pieces = set()
    for i in range(len(protein) - 2):
        pieces.add(protein[i:i + 3])
    return pieces

# Step 3: the more pieces two animals share, the closer they are
def get_distance(protein1, protein2):
    a = get_pieces(protein1)
    b = get_pieces(protein2)
    shared = len(a & b)
    total = len(a | b)
    return 1 - shared / total      # 0 = same, 1 = nothing shared

# Step 4: make the distance table
names = list(sequences)
table = []
for i in range(len(names)):
    row = []
    for j in range(i + 1):
        row.append(get_distance(sequences[names[i]], sequences[names[j]]))
    table.append(row)

distances = DistanceMatrix(names, table)
print(distances)

# Step 5: build and draw the tree
tree = DistanceTreeConstructor().nj(distances)

if "Mus musculus" in sequences:
    tree.root_with_outgroup("Mus musculus")   # mouse is the far-away animal

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax)
plt.title("Cytochrome b tree")
plt.savefig("tree.png", dpi=150)

# Save the sequences and table too
with open("sequences.fasta", "w") as f:
    for name in sequences:
        f.write(">" + name + "\n" + sequences[name] + "\n")

with open("distances.txt", "w") as f:
    f.write(str(distances))

print("Done! Files saved: tree.png, sequences.fasta, distances.txt")
