# Modeling music variations with multiple sequence alignment and profile HMMs
Python routines for multiple sequence alignment, generating profile HMMs, focused on music sequences.

These are my standard set of routines for data-driven analysis of unidimensional music sequences. They contain functions for pairwise alignment (based on BioPython), multiple sequence alignment (using MAFFT or progressive alignment with T-COFFEE), substitution matrix creation from alignments (based on BioPython), consensus generation from an MSA (majority voting, or data fusion), profile HMM modelling of an MSA using Krogh's original model (1994) and HMMER.

##### Example usage

We start with a simple example of protein sequences. First, we create two random 20-symbol sequences from the protein alphabet:

```python
seqA = "".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])
seqB = "".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])

```

Globally pairwise align them and pretty-print them on the terminal:

```python
score,alignment_info=al.pairwiseBetweenSequences(seqA, seqB, match=1, mismatch=-1, gapopen=-0.8, gapext=-0.5, type_="global", returnAlignment=True)

alignment = [alignment_info[0], alignment_info[1]]
al.printMSA(alignment)
print "Alignment score:",score

```

Now, let's add some annotations on the original sequences and visualize them on the terminal:

```python
LABELS = list("-12")
maskA = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
maskB = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
al.printMSAwithMask([seqA,seqB], [maskA,maskB])
```

We now focus on melodic sequences. We willl generate a multiple sequence alignment using MAFFT from a clique of melodic variations (melodies of the same tune family):

```python
sequences, ids = io.readFASTAandIDs("NotAligned/NLBproperSmall/Daar_ging_een_heer_1.fasta")

#ensure the sequences are of the same length first by adding gaps in the end
al.printMSA(al.makeSameLength(sequences)) 

MSA=al.runMafftAlignmentWithSettings(sequences,2,1,method="globalpair",allowshift=False)
print "Aligned sequences:"
al.printMSA(MSA)

```

##### Citation
If you use any part of the code, please cite the following publication:

Bountouridis, D., Brown, D.G., Wiering, F. and Veltkamp, R.C.	Melodic Similarity and Applications Using Biologically-Inspired Techniques. Appl. Sci. 2017, 7, 1242.

