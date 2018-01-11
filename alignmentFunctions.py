from __future__ import division
import os
from os import listdir
from os.path import isfile, isdir, join
import time
import re
import json
import pylev
import numpy as np
from random import shuffle
import matplotlib.pyplot as plt
import csv
import math
import string
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
import random
from types import *
from scipy.interpolate import interp1d
import collections
from termcolor import colored
import matrixFunctions as matrx
import inputOutput as io
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap



BIOALPHABET=list("ARNDCQEGHILKMFPSTWYVBZX") #Protein alphabet



#Generate a substitution matrix for pairwise alignment using an
#array of groups of aligned sequences
def createMatrixFromGroupsOfSequences(groups,name="mymatrix"):
	#maximum number of sequences
	maxG=np.max([len(group) for group in groups])

	#create holder
	sequences=list("*"*maxG)

	#replace gaps with neutral "*" signs
	for group in groups:
		for i,sequence in enumerate(group): sequences[i]+=sequence.replace("-","*")
		
	#export to FASTA and compute LOM
	io.exportGroupOfSequencesToFASTA(sequences,"temp/full1.fasta")
	LOM=matrx.createScoreMatrixFromAlignment("temp/full1.fasta","mylearned-"+name+".matrix")
	
	#normalize the LOM cells
	X=[LOM[tuple] for tuple in LOM]
	mX=np.mean(X)
	stdX=np.std(X)

	#convert LOM to a matrix dictionary and a T-COFFEE matrix
	matrix={}
	MATRIX=np.zeros((23,23))-9
	for tuple in LOM:
		matrix.update({tuple:LOM[tuple] })
		MATRIX[A.index(tuple[0]),A.index(tuple[1])]=((LOM[tuple]-mX)/stdX)*9

	#export
	pickle.dump(matrix,open("mylearned-"+name+".pickle","wb" ) )
	matrx.tcoffeeMakeMatrix(MATRIX,A,"mylearned-"+name+".matrix",-9)
	
	
#Resampling with interpolation a vector
def interpolateResample(vector,size):
	x = np.linspace(0, len(vector), num=len(vector), endpoint=True)
	y = vector
	f = interp1d(x, y,kind='nearest')
	xnew = np.linspace(0, len(vector), num=size, endpoint=True)
	return f(xnew)

#Add gaps to the end of sequences in an MSA so that they are of the same length
def makeSameLength(sequences):
	#maximum sequence length in 
	maxl=np.max([len(sequences[i]) for i in range(0,len(sequences)) ])

	#add gaps
	S=["".join(sequences[i])+"-"*(maxl-len(sequences[i])) for i in range(0,len(sequences))]
	return S

#Naive key estimation of a chord sequence
def findKey(sequence,knownkey=False):
	
	#table of possible chords (function-less) per key
	map=[['A','B','C#','D','E','F#','G#'],
			['A#','C','D','D#','F','G','A'],
			['B','C#','D#','E','F#','G#','A#'],
			['C','D','E','F','G','A','B'],
			['C#','D#','F','F#','G#','A#','C'],
			['D','E','F#','G','A','B','C#'],
			['D#','F','G','G#','A#','C','D'],
			['E','F#','G#','A','B','C#','D#'],
			['F','G','A','A#','C','D','E'],
			['F#','G#','A#','B','C#','D#','F'],
			['G','A','B','C','D','E','F#'],
			['G#','A#','C','C#','D#','F','G']]
	keys=['A','A#','B','C','C#','D','D#','E','F','F#','G','G#' ]
	function=['','m','m','','','m','dim']
	
	#create naive table from functions and chords
	naive_table=[]
	for group in map: 
		t=[v+function[i] for i,v in enumerate(group)]
		naive_table.append(t)

	#Create count histogram
	counts=np.zeros(len(keys))
	for chord in sequence:
		for i,row in enumerate(naive_table):
			if chord in row : counts[i]=counts[i]+1

	#the bin index of the most frequenct bin
	top=np.argmax(counts)

	#check for conflicts
	howmany=0
	for v in counts:
		if v==counts[top] : howmany=howmany+1
	if howmany>=2 or knownkey:
		print "Conflict or known key..."
		newSequence=False
		keyFound=False
		return False
	else: #if no conflicts, return key
		keyFound=keys[top]
		print "Found key:",keyFound
		if knownkey:
			print "Key known"
			keyFound=knownkey
	return keyFound

#compute the distance between two MSAs using the METAL framework
def MSAdistance(msa1,msa2):
	tempfilein1="temp/dis_msa1.fasta"
	tempfilein2="temp/dis_msa2.fasta"
	io.exportGroupOfSequencesToFASTA(msa1,tempfilein1)
	io.exportGroupOfSequencesToFASTA(msa2,tempfilein2)
	command="./metal "+tempfilein1+" "+tempfilein2+" >temp/metal_dis.txt"
	os.system(command)
	txt=readTxtFile("temp/metal_dis.txt")
	sp=txt.split("=")
	return float(sp[1])

#compute the SOP similarity between two MSAs using T-COFFEE
def MSAsimilaritySOP(msa1,msa2):
	tempfilein1="temp/dis_msa1.fasta"
	tempfilein2="temp/dis_msa2.fasta"
	io.exportGroupOfSequencesToFASTA(msa1,tempfilein1)
	io.exportGroupOfSequencesToFASTA(msa2,tempfilein2)
	command="./t_coffee -other_pg aln_compare -al1 "+tempfilein1+" -al2 "+tempfilein2+"  >temp/tcoffee_dis.txt"
	os.system(command)
	txt=readTxtFile("temp/tcoffee_dis.txt")
	txt=txt.split("\n")
	line=txt[2].split(" ")
	values=[]
	for item in line:
		try:
			values.append(float(item))
		except ValueError:
			i=1
	score=( values[-1] )/100
	return score

#perform inner metric analysis on an vector of onset values
def innerMetricAnalysis(onset_vector,type="metrical"):
	np.savetxt('temp/ima.in', np.array(onset_vector), delimiter=',',fmt='%d',) 
	command="./IMA -m metrical -f 'temp/ima.in' -o 'temp/ima.out' "
	os.system(command)
	
	if type=="spectral":
		command="./IMA -m "+type+" -f 'temp/ima.in' -o 'temp/ima_s.out' "
		os.system(command)	
	
	ima=np.loadtxt('temp/ima.out')
	if type=="spectral":
		imas=np.loadtxt('temp/ima_s.out')
		imas=imas/np.max(imas)
		new_imas = np.floor(np.delete(imas, np.where(ima==0))/0.044)
		return new_imas
	else:
		ima=ima/np.max(ima)
		new_ima = np.floor(np.delete(ima, np.where(ima==0))/0.044)
		return new_ima

#compute the naive concensus sequence from an alignment of sequences
def consensusFromSequencesCnf(sequences,threshold,ambiguous,require_multiple):
	io.exportGroupOfSequencesToFASTA(sequences,"temp/forconsensus.fasta")
	c_align = AlignIO.read("temp/forconsensus.fasta", "fasta")
	summary_align = AlignInfo.SummaryInfo(c_align)
	concensus=summary_align.dumb_consensus(threshold=threshold, ambiguous=ambiguous, consensus_alpha=None, require_multiple=require_multiple)
	return concensus

#run MAFFT alignment on a group of sequences without a substitution matrix
def runMafftAlignmentWithSettings(inputSequences,gap_open,gap_ext,method="globalpair",allowshift=False):
	tempfile="temp/tempseq.fasta"
	tempfileout="temp/tempmsa.fasta"
	tempfileout2="temp/tempmsa2.fasta"
	io.exportGroupOfSequencesToFASTA(inputSequences,tempfile)
	if allowshift==True:
		command="mafft-mac/mafft.bat --quiet --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000 --allowshift  --text "+tempfile+"		> "+tempfileout2
	else:
		command="mafft-mac/mafft.bat --quiet --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000   --text "+tempfile+"		> "+tempfileout2
	os.system(command)
	output=io.readFASTA(tempfileout2)
	return output

#run MAFFT alignment on a group of sequences with a substitution matrix
def runMafftAlignmentWithMatrix(inputSequences,gap_open,gap_ext,method="globalpair",matrixfile=""):
	tempfile="temp/tempseq.fasta"
	tempfileout="temp/tempmsa.fasta"
	tempfileout2="temp/tempmsa2.fasta"
	io.exportGroupOfSequencesToFASTA(inputSequences,tempfile)
	command="mafft-mac/mafft.bat --quiet  --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000 --aamatrix "+matrixfile+" "+tempfile+"		> "+tempfileout2
	os.system(command)
	output=io.readFASTA(tempfileout2)
	return output

#compute the consensus sequence from an MSA using the data fusion technique
def fusionConsensus(sequences,tempin,tempout):
	exportGroupOfSequencesToString(sequences,tempin)

	#the domain needs to be established first
	Alphabet=list(string.ascii_uppercase)
	for k in range(0,len(sequences[0])):
		if sequences[0][k] in Alphabet:
			domainfile="abc.dom"
			break
		if sequences[0][k] in [0,1] or sequences[0][0] in ["0","1"]:
			domainfile="binary.dom"
			break
	command='./fusion 5 "'+tempin+'"  "'+domainfile+'" "'+tempout+'" "'+tempout+'.acc"  >"temp/fusionverbose.txt"'
	os.system(command)
	consensus=readTxtFile(tempout)
	return consensus

#pairwise alignment between two sequences with input gap settings and substitution matrix
def pairwiseBetweenSequencesWithMatrix(seqA, seqB, gapopen=-0.8, gapext=-0.5, matrix=[], type_="global",returnAlignment=False):
    maxscore = -1000000
    if type_ == "global":
    	if not returnAlignment:
        	return pairwise2.align.globalds(seqA, seqB, matrix, gapopen, gapext, score_only=True),False
        else:
	        for a in pairwise2.align.globalds(seqA, seqB, matrix, gapopen, gapext, score_only=False):
	            score = a[2]
	            if score > maxscore:
	                winSeq = seqB
	                maxscore = score
	                winAli = a
	            break
	        return score,winAli
    if type_ == "local":
    	if not returnAlignment:
        	return pairwise2.align.localdx(seqA, seqB, matrix, gapopen, gapext, score_only=True),False
        else:
	        for a in pairwise2.align.localdx(seqA, seqB, matrix, gapopen, gapext, score_only=False):
	            score = a[2]
	            if score > maxscore:
	                winSeq = seqB
	                maxscore = score
	                winAli = a
	            break
	        return score,winAli

#pairwise alignment between two sequences with input gap settings match mismatch values  
def pairwiseBetweenSequences(seqA,seqB, match=1,mismatch=-1,gapopen=-0.8,gapext=-0.5,type_="global",returnAlignment=False):
	maxscore=-1000000
	if type_=="global":
		if not returnAlignment:
			return pairwise2.align.globalms(seqA, seqB, match,mismatch,gapopen,gapext,score_only=True),False
		else:
			for a in pairwise2.align.globalms(seqA, seqB, match,mismatch,gapopen,gapext):
				score=a[2]
				if score>maxscore:
					winSeq=seqB
					maxscore=score
					winAli=a
				break
			return score,winAli
	if type_=="local":
		if not returnAlignment:
			return pairwise2.align.localms(seqA, seqB, match,mismatch,gapopen,gapext,score_only=True),False
		else:
			for a in pairwise2.align.localms(seqA, seqB, match,mismatch,gapopen,gapext):
				score=a[2]
				if score>maxscore:
					winSeq=seqB
					maxscore=score
					winAli=a
				break
			return score,winAli

#Print a multiple sequence alignment with mask (corresponding to patten labels 0,1,2) on the terminal
def printMSAwithMask(MSA,mask):
    Alphabet=list("123456789")

    #create the palette of possible foreground background combinations
    col1=[("red",[]),("green",[]),("yellow",[]),("blue",[]),("magenta",[]),("cyan",[]),("white",[])]
    colors=col1
    for i in range(0,len(col1)):
        color1=col1[i][0]
        for j in range(0,len(col1)):
            color2=col1[j][0]
            if color1!=color2:  colors.append((color1,"on_"+color2))
            if len(colors)>=30: break

    #assign color to symbols
    for i,sequence in enumerate(MSA):
        text=""
        for j,c in enumerate(sequence):
            #print Alphabet.index(c)
            if mask[i][j] not in Alphabet:
                CL=(("grey",[]))
            else:
                CL=colors[Alphabet.index(mask[i][j])]
            if len(CL[1])<1:text+=colored(c,CL[0])
            else:text+=colored(c,CL[0],CL[1])
        print text

#Print an MSA mask on the terminal
def printMSAmask(MSA):
    Alphabet=["-"] + list("123")+list("ABCDEFG") +list("abcdefg")+list("HIJKLMN") +list("hijklmn")

    #create the palette of possible foreground background combinations
    colors=[("grey",[]),("grey","on_red"),("grey","on_green"),("grey","on_yellow"),("grey","on_blue"),("grey","on_magenta"),("grey","on_cyan"),("white",[])]
    col1=[("grey",[]),("red",[]),("green",[]),("blue",[]),("yellow",[]),("magenta",[]),("cyan",[]),("white",[])]
    colors=col1
    for i in range(0,len(col1)):
        color1=col1[i][0]
        for j in range(0,len(col1)):
            color2=col1[j][0]
            if color1!=color2:  colors.append((color1,"on_"+color2))
            if len(colors)>=42: break

    #assign colors to symbols
    for sequence in MSA:
        text=""
        for c in sequence :
            CL=colors[Alphabet.index(c)]
            if len(CL[1])<1:text+=colored(c,CL[0])
            else:text+=colored(c,CL[0],CL[1])
        print text

#Print an MSA on the terminal
def printMSA(MSA):
    Alphabet=list("-ARNDCQEGHILKMFPSTWYVBZX12345678*")

    #create the palette of possible foreground background combinations
    col1=[("grey",[]),("red",[]),("green",[]),("yellow",[]),("blue",[]),("magenta",[]),("cyan",[]),("white",[])]
    colors=col1
    for i in range(0,len(col1)):
        color1=col1[i][0]
        for j in range(0,len(col1)):
            color2=col1[j][0]
            if color1!=color2:  colors.append((color1,"on_"+color2))
            if len(colors)>=32: break

    #assign colors to symbols
    for sequence in MSA:
        text=""
        for c in sequence :
            CL=colors[Alphabet.index(c)]
            if len(CL[1])<1:text+=colored(c,CL[0])
            else:text+=colored(c,CL[0],CL[1])
        print text

#latex format of multiple sequence alignment with mask
def plotOneAlignmentWithMaskLatex(sequences,mask,output):
    colors=["OliveGreen","red","cyan", "blue","Melon","SkyBlue"]
    s="\\pbox{20cm}{"
    for i,seq in enumerate(sequences):
        s+="\\texttt{"
        for j,char in enumerate(seq):
            if mask[i][j]!="-":
                color=colors[int(mask[i][j])-1]
                s+="\\textcolor{"+color+"}"
            s+="{"+char+"}"
        s+="}"
        s+="\\\\ \n"
    s+="}"
    return s
   
#seaborn heatmap of multiple sequence alignment with mask
def plotOneAlignmentWithMask(sequences,mask,output,blackAndWhite=False):
    map1=["-","1","2","3","4","5","6","7","8","9"]
    size=len(sequences)
    length=len(sequences[0])
    M = np.zeros((size,length))
    for i, seq in enumerate(mask):
        for j in range(length):
            M[i][j] = map1.index(seq[j])
    corr = pd.DataFrame(M,  columns=range(length))
    sns.set_context("notebook", font_scale=0.8, rc={"lines.linewidth": 1.0})
    sns.set_style({'font.family': 'serif', 'font.serif': ['Times New Roman']})

    # Set up the matplotlib figure
    _, ax = plt.subplots(figsize=(15,15*size/length))
    cmap=sns.light_palette((210, 90, 60), input="husl")
    flatui = ["#ffffff", "#a2cffe","#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    cmap=sns.cubehelix_palette(dark=0.3, light=1, as_cmap=True,rot=.1)
    sns.heatmap(corr,cbar=False,  ax=ax, annot=True,fmt=".1f",linewidths=.5,cmap=cmap)

    ALPHA=["-"+BIOALPHABET]

    s=np.array(sequences).flatten()
    for i,text in enumerate(ax.texts):
        if i<len(s):
            if s[i]!="-":text.set_text(s[i])
            else: text.set_text(" ")
            text.set_size(8)
            text.set_weight('bold')
    plt.savefig(output,format="pdf")
    plt.close()

#plot alignment with our without stability histogram
def plotAlignmentWithNaiveStability(sequences,output,labels=False):
    Alphabet=BIOALPHABET[:]

    Alphabet=list(set([s for s in np.array(sequences).flatten()]))
    Alphabet.pop(Alphabet.index("-"))
    shuffle(Alphabet)
    print Alphabet


 
    #convert to image
    image=np.zeros((len(sequences),len(sequences[0])))
    x=0
    for sequence in sequences:
        y=0
        for s in sequence:
            if s=="-": image[x][y]=0
            else:
                image[x][y]=Alphabet.index(s)+1
            y+=1
        x+=1
	fig = plt.figure(figsize=(15, 8))
	gs1 = gridspec.GridSpec(4, 4)
	gs1.update( hspace=0.00,wspace=0.00)
	axes0 = plt.subplot(gs1[:-1, :])
	axes1 = plt.subplot(gs1[-1, :])

	colors_ = [(256,256,256),(32,185,171),(194,156,21),(25,147,136),(8,120,115),(0,105,101),(0,60,72),(249,115,6),(252,90,80),(233,187,25)]  # R -> G -> B
  	C=[]
  	for c in colors_+colors_[1:]:
  		r,g,b=c
  		C.append((r/256.0,g/256.0,b/256.0))
  	
  	
    cm = LinearSegmentedColormap.from_list("mycmap", C, N=len(C))


    conf_arr = image
    axes0.set_aspect(1)
    axes0.imshow(image, cmap=cm, interpolation='nearest')
    width = len(conf_arr)
    height = len(conf_arr[0])
    for x in xrange(width):
        for y in xrange(height):
            axes0.annotate(sequences[x][y], xy=(y, x), 
                        horizontalalignment='center',
                        verticalalignment='center',size=10,color="white")
    axes0.axis('off')
    s=np.array(sequences).flatten()
    for i,text in enumerate(axes0.texts):
        if s[i]=="-":text.set_color("black")
       
    # The second axis is the "stability"
	axes1.bar(range(0,len(sequences[0])), mostFrequentSymbolPerColumn(sequences),color='gray',width=1.0,linewidth=0)
	axes1.set_xlim([0,len(sequences[0])])
    #plt.savefig(output,format="pdf")
    axes1.spines['right'].set_visible(False)
    axes1.spines['top'].set_visible(False)
    axes1.spines['bottom'].set_visible(False)

    plt.show()
    plt.close()

#compute the frequency of the most frequence symbol per column in an MSA
def mostFrequentSymbolPerColumn(MSA):
	numOfsequences=len(MSA)
	lengthOfMSA=len(MSA[0])
	a=np.array(MSA)
	
	ratios=[]
	for i in xrange(lengthOfMSA):
		column=a[:,i]
		gapless=column[np.where(column!="-")]
		counts=collections.Counter(column[np.where(column!="-")])
		#percentageOfgaps.append(len(np.where(column=="-")[0])/numOfsequences)
		mx=np.max([counts[char] for char in counts])
		ratio=mx/numOfsequences
		ratios.append(ratio)
	return ratios