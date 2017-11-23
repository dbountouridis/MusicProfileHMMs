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