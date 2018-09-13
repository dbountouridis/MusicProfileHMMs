# Modeling music variations with multiple sequence alignment and profile HMMs
Python routines for data-driven analysis of unidimensional music sequences. Some of the functions include pairwise alignment (based on [BioPython](http://biopython.org/)), multiple sequence alignment (using [MAFFT](https://mafft.cbrc.jp/alignment/software/) or progressive alignment with [T-COFFEE](http://www.tcoffee.org/)), substitution matrix creation from alignments (based on BioPython), consensus generation from an MSA (majority voting), profile HMM modelling of an MSA using Krogh's original architecture (1994) and [HMMER](http://hmmer.org/) (modified for music sequences).

### Example usage

We will go through the _demo.py_ script to get a glimpse of the routines' functionalities. Before dealing with music sequences, we start with a simple example of aligning protein sequences. First, we create two random 20-symbol sequences from the protein alphabet:

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

The output should look something like this:

![alt text](https://github.com/dbountouridis/MusicProfileHMMs/blob/master/images/1a.png "Alignment1")

Add some random annotations on the original sequences and visualize them on the terminal:

```python
LABELS = list("-12")
maskA = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
maskB = "".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
al.printMSAwithMask([seqA,seqB], [maskA,maskB])
```

We now focus on the more interesting stuff, melodic sequences. We willl generate a multiple sequence alignment using MAFFT from a clique of melodic variations (melodies of the same [MTC](http://www.liederenbank.nl/mtc/) tune family). The melodies have already been converted into alphabetic sequences of pitch intervals:

```python
sequences, ids = io.readFASTAandIDs("NotAligned/NLBproperSmall/Daar_ging_een_heer_1.fasta")

MSA = al.runMafftAlignmentWithSettings(sequences, 2, 1, method="globalpair", allowshift=False)
al.printMSA(MSA)

```

The output should look something like this:

![alt text](https://github.com/dbountouridis/MusicProfileHMMs/blob/master/images/4a.png "Alignment3")

The next step is to model the multiple sequene alignment using a profile HMM. But first we need to compute the global frequency of the symbols (pitch intervals) in the MTC dataset:

```python
directory = "NotAligned/NLBproperSmall"
files = io.filesInPath(directory)

allPossibleSymbols = "".join(["".join(["".join(seq) for seq in io.readFASTA(directory+"/"+file)]) for file in files]).replace("-", "")
counts = collections.Counter(allPossibleSymbols)  
alphabet = [symbol for symbol in counts]  
counts = np.array([counts[symbol] for symbol in counts])  

EmissionProbabilities = dict( zip(alphabet, counts/np.sum(counts)))  # convert to dictionary

```

Create a profile-HMM from the MSA:

```python
sw = 10
pct = 1
pec = 1/12.
grms = 0.95
model, columnsIsMatchState = krogh.profileHMM(MSA, alphabet=alphabet,gapRatiotoMatchStates=grms, pseudoCountTransition=pct,
	sequenceWeight=sw, pseudoEmissionCount=pec, 
	plot_=False, uniformInsertionProbs=True, nullEmission=EmissionProbabilities)

# Train the profile HMM on the MSA to re-adjust the HMM parameters
multiplication = 10
di = 0.5
ei = 0.5
s_th = 0.000001
model.train([sequence for sequence in MSA], max_iterations = 500, distribution_inertia = di,edge_inertia = ei, stop_threshold = s_th, algorithm = 'baum-welocalh')

```

Finally, comprare sequences of a random family plus a sequence from the original MSA to the profile HMM:

```python
familyIndex = int(random.random()*len(files))
querySequences = ["".join(seq).replace("-", "") for seq in io.readFASTA(directory+"/"+files[familyIndex])]+["".join(MSA[0])]
queryIds = [files[familyIndex]] * (len(querySequences)-1)+["Daar_ging_een_heer_1.fasta"]

scores = [krogh.compareSequenceToProfHMM(model, query)[0] for query in querySequences]
labels = ["Daar_ging_een_heer_1.fasta" == filename for filename in queryIds]
print "Average precision:", average_precision_score(labels, scores)
```

For modelling music variations using the HMMER architecture of profile-HMMs please email me at d.bountouridis@gmail.com. HMMER is generally faster but its accuracy is not much better than the Krogh architecture (for music sequences).


### Citation
If you use any part of the code, please cite the following publication:

Bountouridis, D., Brown, D.G., Wiering, F. and Veltkamp, R.C.	_Melodic Similarity and Applications Using Biologically-Inspired Techniques_. Appl. Sci. 2017, 7, 1242.

or 

Bountouridis D. _Music Information Retrieval Using Biologically Inspired Techniques
_ PhD Thesis, Utrecht University, 2018.

