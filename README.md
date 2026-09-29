# Primate Phylogeny from Cytochrome b

A small Python project I made while learning bioinformatics.

## Question
Can protein sequences show how closely related some mammals are?

## Data
Cytochrome b protein sequences (about 380 amino acids each) from NCBI for 7 animals: human, chimpanzee, gorilla, orangutan, macaque, mouse and cow.

## Method
1. Downloaded the sequences with Biopython (Entrez)
2. Cut each protein into pieces of 3 letters
3. Measured how many pieces two animals share (distance = 1 - shared / total)
4. Built a distance table
5. Made a tree with the Neighbor-Joining method and drew it

## Results
- Human and chimpanzee were the closest pair (distance 0.27)
- Human and gorilla came next (0.28)
- Mouse and cow were far from the primates

## Limits
- This method does not line up the sequences, so it is simple but not very exact
- Macaque came out further from human than expected
- Cow came out closer to mouse than to primates, which is not right biologically
- Only one gene and 7 animals were used

## What I learned
- How to download sequences from NCBI with Biopython
- How a distance table becomes a tree
- Why simple methods can give odd results

## Next steps
- Use multiple sequence alignment (Clustal Omega or MUSCLE)
- Use more genes and more animals
- Try another tree method and compare

## Tools
Python, Biopython, Matplotlib, Google Colab.
