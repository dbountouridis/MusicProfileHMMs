# Modeling music variations with multiple sequence alignment and profile HMMs
Python routines for data-driven analysis of unidimensional music sequences. Some of the functions include pairwise alignment (based on [BioPython](http://biopython.org/)), multiple sequence alignment (using [MAFFT](https://mafft.cbrc.jp/alignment/software/) or progressive alignment with [T-COFFEE](http://www.tcoffee.org/)), substitution matrix creation from alignments (based on BioPython), consensus generation from an MSA (majority voting), profile HMM modelling of an MSA using Krogh's original architecture (1994) and [HMMER](http://hmmer.org/) (modified for music sequences).

##### Example usage

Before dealing with music sequences, we start with a simple example of aligning protein sequences. First, we create two random 20-symbol sequences from the protein alphabet:

```python
BIOALPHABET = list("ARNDCQEGHILKMFPSTWYVBZX")

seqA = "".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])
seqB = "".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])

```

Globally pairwise align them and pretty-print them on the terminal:

```python
score, alignment_info = al.pairwiseBetweenSequences(seqA, seqB, match=1, mismatch=-1, gapopen=-0.8, gapext=-0.5, type_="global", returnAlignment=True)

alignment = [alignment_info[0], alignment_info[1]]
al.printMSA(alignment)

```

Add some annotations on the original sequences and visualize them on the terminal:

```python
LABELS = list("-12")
maskA = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
maskB = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
al.printMSAwithMask([seqA,seqB], [maskA,maskB])
```

We now focus on melodic sequences. We willl generate a multiple sequence alignment using MAFFT from a clique of melodic variations (melodies of the same [MTC](http://www.liederenbank.nl/mtc/) tune family). The melodies have already been converted into alphabetic sequences of pitch intervals:

```python
sequences, ids = io.readFASTAandIDs("NotAligned/NLBproperSmall/Daar_ging_een_heer_1.fasta")

MSA = al.runMafftAlignmentWithSettings(sequences, 2, 1, method="globalpair", allowshift=False)
al.printMSA(MSA)

```

##### Citation
If you use any part of the code, please cite the following publication:

Bountouridis, D., Brown, D.G., Wiering, F. and Veltkamp, R.C.	_Melodic Similarity and Applications Using Biologically-Inspired Techniques_. Appl. Sci. 2017, 7, 1242.

