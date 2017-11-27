from __future__ import division
import os
import time
import math
import numpy as np
from random import shuffle
import string
import random
from types import *
import collections
import alignmentFunctions as al
import inputOutput as io
import kroghProfileHMM as krogh

BIOALPHABET=list("ARNDCQEGHILKMFPSTWYVBZX") #Protein alphabet


#create two random sequences from the protein alphabet
seqA="".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])
seqB="".join([BIOALPHABET[int(random.random()*len(BIOALPHABET))] for i in range(20)])

#globally align them
score,alignment_info=al.pairwiseBetweenSequences(seqA,seqB, match=1,mismatch=-1,gapopen=-0.8,gapext=-0.5,type_="global",returnAlignment=True)

#Beautifully showing the alignment on the terminal
alignment=[alignment_info[0],alignment_info[1]]
al.printMSA(alignment)
print "Alignment score:",score

#let's add some annotations on the original sequences and visualize them on the terminal
LABELS=list("-12")
maskA="".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
maskB="".join([LABELS[int(random.random()*len(LABELS))] for i in range(20)])
al.printMSAwithMask([seqA,seqB],[maskA,maskB])

#let's do some multiple sequence alignment
#first let's load sequences corresponding to melodies of the same tune family 
sequences,ids=io.readFASTAandIDs("NotAligned/NLBproperSmall/Daar_ging_een_heer_1.fasta")

print "Unaligned sequences:"
al.printMSA(al.makeSameLength(sequences)) #ensure the sequences are of the same length first by adding gaps in the end

#align them using MAFFT
MSA=al.runMafftAlignmentWithSettings(sequences,2,1,method="globalpair",allowshift=False)
print "Aligned sequences:"
al.printMSA(MSA)

#Read MIDI file
Notes,Durations,Onsets,minduration,PitchIntervals=io.readMIDIfile("AM-13A-2.mid")

#Create a Krogh profile HMM from a folk tune family
#First we need the emission probabilities gathered from the whole dataset
directory="NotAligned/NLBproperSmall"
files=io.filesInPath(directory)
allPossibleSymbols="".join(["".join(["".join(seq) for seq in io.readFASTA(directory+"/"+file)]) for file in files]).replace("-","")
counts=collections.Counter(allPossibleSymbols) #count
alphabet=[symbol for symbol in counts] #all possible symbols
counts=np.array([counts[symbol] for symbol in counts]) #corresponding counts
EmissionProbabilities = dict(zip(alphabet, counts/np.sum(counts) )) #convert to dictionary

#Bake the Krogh profile HMM model
sw=10
pct=1
pec=1/12.
grms=0.95
model,columnsIsMatchState=krogh.profileHMM(MSA,alphabet=alphabet,gapRatiotoMatchStates=grms,pseudoCountTransition=pct,sequenceWeight=sw,pseudoEmissionCount=pec,plot_=False,uniformInsertionProbs=True,nullEmission=EmissionProbabilities)

#Train it on the MSA to readjust the HMM parameters
multiplication=10
di=0.5
ei=0.5
s_th=0.000001
model.train([sequence for sequence in MSA],max_iterations=500,distribution_inertia=di, edge_inertia=ei,stop_threshold=s_th,algorithm='baum-welocalh')

#Compare sequences to the profile HMM


	






	