""" A collection of functions for alignment-related tasks. 

The collection includes generating substitution matrices from alignments,
computing the distance between MSAs, computing MSAs (using MAFFT, T-COFFEE) and so on.

Example:
	See example.py

Todo:
	* Convert these functions into an alignment class

"""

from __future__ import division

import collections
import csv
import inputOutput as io
import json
import math
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import pylev
import random
import re
import string
import time

from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC
from Bio.pairwise2 import format_alignment
from matplotlib.colors import LinearSegmentedColormap
from os import listdir
from os.path import isdir
from os.path import isfile
from os.path import join
from random import shuffle
from scipy.interpolate import interp1d
from termcolor import colored
from types import *

__author__ = 'Dimitrios  Bountouridis'

# Protein alphabet
BIOALPHABET = list("ARNDCQEGHILKMFPSTWYVBZX") 

def createMatrixFromGroupsOfSequences(groups, name="mymatrix"):
	""" Generate a substitution matrix for pairwise alignment using an
	array of groups (cliques) of alignments.

	This function is used in the case where a matrix needs to be computed
	from a number of cliques of aligned sequences. For example, the aligned
	Liederenbank dataset. The function generates a matrix as a dictionary
	e.g. {(A,B):float,(A,C):float}} and so on (similar to the biopython pairwise matrix).
	It also generates a matrix in the format for T-COFFEE (which needs ints).

	Args:
		groups (array): an array of arrays containing aligned strings
		name (string): the output matrix file name


	"""

	# Maximum number of sequences
	maxG = np.max([len(group) for group in groups])

	# Create holder
	sequences = list("*"*maxG)

	# Replace gaps with neutral "*" signs
	for group in groups:
		for i,sequence in enumerate(group): sequences[i]+=sequence.replace("-","*")
		
	# Export to FASTA and compute LOM
	io.exportGroupOfSequencesToFASTA(sequences,"temp/full1.fasta")
	LOM = createScoreMatrixFromAlignment("temp/full1.fasta","mylearned-"+name+".matrix")
	
	# Normalize LOM 
	X = [LOM[tuple] for tuple in LOM]
	mX = np.mean(X)
	stdX = np.std(X)

	# Convert LOM to a matrix dictionary and a T-COFFEE matrix
	matrix = {}
	MATRIX = np.zeros((lent(BIOALPHABET),len(BIOALPHABET)))-9
	for tuple in LOM:
		matrix.update({tuple:LOM[tuple] })
		MATRIX[BIOALPHABET.index(tuple[0]), BIOALPHABET.index(tuple[1])] = ((LOM[tuple]-mX)/stdX)*9

	# Export
	pickle.dump(matrix,open("mylearned-"+name+".pickle","wb" ) )
	tcoffeeMakeMatrix(MATRIX,"mylearned-"+name+".matrix",-9)
	
def makeSameLength(sequences):
	""" Add gaps to the end of sequences so that they are of the same length.

	Args:
		sequences (array of strings): strings to be adjusted for length

	Returns:
		S: sequences of same length


	"""

	maxl = np.max([len(sequences[i]) for i in range(0,len(sequences)) ])
	S = ["".join(sequences[i])+"-"*(maxl-len(sequences[i])) for i in range(0,len(sequences))]
	return S


def sumOfPairsDistanceMSA(msa1, msa2):
	""" Compute the SOP distance between two MSAs using T-COFFEE

	The MSAs should correspond to alignments of the same sequences.

	Args:
		msa1 (array): array of aligned strings
		msa1 (array): array of aligned strings

	Returns:
		score (float): distance between two MSAs


	"""

	tempfilein1 = "temp/dis_msa1.fasta"
	tempfilein2 = "temp/dis_msa2.fasta"
	io.exportGroupOfSequencesToFASTA(msa1,tempfilein1)
	io.exportGroupOfSequencesToFASTA(msa2,tempfilein2)
	command = "./t_coffee -other_pg aln_compare -al1 "+tempfilein1+" -al2 "+tempfilein2+"  >temp/tcoffee_dis.txt"
	os.system(command)
	txt = readTxtFile("temp/tcoffee_dis.txt")
	txt = txt.split("\n")
	line = txt[2].split(" ")
	values = []
	for item in line:
		try:
			values.append(float(item))
		except ValueError:
			i = 1
	score = values[-1]
	return score

def consensusFromSequencesCnf(sequences, threshold, ambiguous, require_multiple):
	""" Compute the naive concensus sequence from an alignment of sequences.

	Args:
		sequences (array): the array of aligned sequences
		threshold (float): the value of column agreement under which a single
						letter will be asigned to the column
		ambiguous (char): a char letter that will be assigned to a column
						for an agreement lower than the threshold
		require_multiple: see biopython.org/DIST/docs/api/Bio.Align.AlignInfo.SummaryInfo-class.html

	Return:
		consensus (str): the consensus sequence


	"""

	io.exportGroupOfSequencesToFASTA(sequences,"temp/forconsensus.fasta")
	c_align = AlignIO.read("temp/forconsensus.fasta", "fasta")
	summary_align = AlignInfo.SummaryInfo(c_align)
	concensus = summary_align.dumb_consensus(threshold = threshold, 
		ambiguous = ambiguous, 
		consensus_alpha = None, 
		require_multiple = require_multiple)
	return concensus

def runMafftAlignmentWithSettings(inputSequences, gap_open, gap_ext, method = "globalpair", allowshift = False):
	""" Run MAFFT alignment on a group of sequences without a substitution matrix.

	This function requires a proper path to the mafft binary.

	Args:
		inputSequences (array): an array of strings of sequences to be aligned
		gap_open (float): gap open penalty
		gap_ext (float): gap extent penalty
		method (str): MAFFT's alignment strategy
		allowshift (bool):  see MAFFT documentation

	Returns:
		output (array): aligned strings of sequences

	"""


	tempfile = "temp/tempseq.fasta"
	tempfileout = "temp/tempmsa.fasta"
	tempfileout2 = "temp/tempmsa2.fasta"
	io.exportGroupOfSequencesToFASTA(inputSequences,tempfile)
	
	if allowshift==True:
		command = "mafft-mac/mafft.bat --quiet --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000 --allowshift  --text "+tempfile+" > "+tempfileout2
	else:
		command = "mafft-mac/mafft.bat --quiet --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000   --text "+tempfile+" > "+tempfileout2
	os.system(command)
	
	output = io.readFASTA(tempfileout2)
	return output

def runMafftAlignmentWithMatrix(inputSequences, gap_open, gap_ext, method="globalpair", matrixfile=""):
	""" Run MAFFT alignment on a group of sequences with a substitution matrix.

	This function requires a proper path to the mafft binary.

	Args:
		inputSequences (array): an array of strings of sequences to be aligned
		gap_open (float): gap open penalty
		gap_ext (float): gap extent penalty
		method (str): MAFFT's alignment strategy
		matrixfile (str): path to an appropriate matrix file for MAFFT

	Returns:
		output (array): aligned strings of sequences

	"""
	tempfile = "temp/tempseq.fasta"
	tempfileout = "temp/tempmsa.fasta"
	tempfileout2 = "temp/tempmsa2.fasta"

	io.exportGroupOfSequencesToFASTA(inputSequences,tempfile)
	command = "mafft-mac/mafft.bat --quiet  --op "+str(gap_open)+" --ep "+str(gap_ext)+" --"+method+" --maxiterate 1000 --aamatrix "+matrixfile+" "+tempfile+"		> "+tempfileout2
	os.system(command)

	output=io.readFASTA(tempfileout2)
	return output


def pairwiseBetweenSequencesWithMatrix(seqA, seqB, gapopen=-0.8, gapext=-0.5, matrix=[], type_="global", returnAlignment=False):
	""" Pairwise alignment between two sequences with input gap settings and substitution matrix.
	
	This function can return the score of the alignment and the alignment itself. But the alignment itself
	is computationaly expensive, so it is suggested to skip  it when only the score of the alignment is of 
	importance.

	Args:
		seqA (str): sequence to be aligned
		seqB (str): sequence to be aligned
		gapopen (float): gap open penatly
		gapext (float): gap extend penalty
		matrix (dict): a matrix dictionary of the form {(A,B):float}
		type_ (str): type of alignment (global or local)
		returnAlignment (bool): whether to return the alignment itself or just the score

	Returns:
		param1 (float): the alignment score
		param2 (): either the alignment or False

	"""
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

 
def pairwiseBetweenSequences(seqA, seqB, match=1, mismatch=-1, gapopen=-0.8, gapext=-0.5, type_="global", returnAlignment=False):
	""" Pairwise alignment between two sequences with input gap settings and match/mismatch values.
	
	This function can return the score of the alignment and the alignment itself. But the alignment itself
	is computationaly expensive, so it is suggested to skip  it when only the score of the alignment is of 
	importance.

	Args:
		seqA (str): sequence to be aligned
		seqB (str): sequence to be aligned
		match (float): the char match score
		mismatch (float): the char mismatch score
		gapopen (float): gap open penatly
		gapext (float): gap extend penalty
		type_ (str): type of alignment (global or local)
		returnAlignment (bool): whether to return the alignment itself or just the score

	Returns:
		param1 (float): the alignment score
		param2 (): either the alignment or False

	"""
	maxscore = -1000000
	if type_=="global":
		if not returnAlignment:
			return pairwise2.align.globalms(seqA, seqB, match,mismatch,gapopen,gapext,score_only=True),False
		else:
			for a in pairwise2.align.globalms(seqA, seqB, match,mismatch,gapopen,gapext):
				score = a[2]
				if score>maxscore:
					winSeq = seqB
					maxscore = score
					winAli = a
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


def printMSAwithMask(MSA,mask):
	""" A pretty terminal print of an MSA with a corresponding mask (labels for each char).

	This function is typically used for investigating whethere an MSA results to an 
	alignment of patterns (labels).

	Example:
		MSA: AA-BC-AF   mask: 11-12-34
			 AC-BB--F         12-22--4

	Args:
		MSA (array): an array of aligned string sequences
		mask (array): an array of aligned string labels
	
	
	"""

	Alphabet = list("123456789")

	# Create the palette of possible foreground background combinations
	col1 = [("red",[]),("green",[]),("yellow",[]),("blue",[]),("magenta",[]),("cyan",[]),("white",[])]
	colors = col1
	for i in range(0,len(col1)):
		color1 = col1[i][0]
		for j in range(0,len(col1)):
			color2 = col1[j][0]
			if color1!=color2:  colors.append((color1,"on_"+color2))
			if len(colors)>=30: break

	# Assign color to symbols
	for i,sequence in enumerate(MSA):
		text = ""
		for j,c in enumerate(sequence):
			if mask[i][j] not in Alphabet:
				CL = (("grey",[]))
			else:
				CL = colors[Alphabet.index(mask[i][j])]
			if len(CL[1])<1:text+=colored(c,CL[0])
			else: text+=colored(c,CL[0],CL[1])
		print text


def printMSA(MSA):
	""" A pretty print of an MSA on the terminal.

	Args:
		MSA: (array): an array of aligned string sequences
	
	"""

	Alphabet = list("-ARNDCQEGHILKMFPSTWYVBZX12345678*")

	# Create the palette of possible foreground background combinations
	col1 = [("grey",[]),("red",[]),("green",[]),("yellow",[]),("blue",[]),("magenta",[]),("cyan",[]),("white",[])]
	colors = col1
	for i in range(0,len(col1)):
		color1 = col1[i][0]
		for j in range(0,len(col1)):
			color2 = col1[j][0]
			if color1!=color2:  colors.append((color1,"on_"+color2))
			if len(colors)>=32: break

	for sequence in MSA:
		text = ""
		for c in sequence :
			CL = colors[Alphabet.index(c)]
			if len(CL[1])<1:text+=colored(c,CL[0])
			else:text+=colored(c,CL[0],CL[1])
		print text


def tcoffeeAlignment(sequences,go,ge):
	""" The standard T-COFFEE alignment.

	For more info: http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
	
	Args:
		sequences (array): array of strings to be aligned
		go (int): gap open penalty
		ge (int): gap extent penalty

	Returns:
		output (array): aligned sequences


	"""

	io.exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command = "./t_coffee 'temp/tcoffee-in.fasta' -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output = io.readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeAlignmentWithID(sequences,go,ge):
	""" The standard T-COFFEE alignment with ID matrix (ones in the diagonal)

	For more info: http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
	
	Args:
		sequences (array): array of strings to be aligned
		go (int): gap open penalty
		ge (int): gap extent penalty

	Returns:
		output (array): aligned sequences


	"""

	io.exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command = "./t_coffee 'temp/tcoffee-in.fasta' -matrix ID  -gapopen = "+str(go)+" -gapext = "+str(ge)+" -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output = io.readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeAlignmentWithMatrix(sequences,go,ge,matrixfile):
	""" The standard T-COFFEE alignment with custom matrix 

	For more info: http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
	
	Args:
		sequences (array): array of strings to be aligned
		go (int): gap open penalty
		ge (int): gap extent penalty
		matrixfile (str): file path to an appropriate matrix file

	Returns:
		output (array): aligned sequences

	"""

	io.exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command = "./t_coffee 'temp/tcoffee-in.fasta' -in = X"+matrixfile+"  -gapopen = "+str(go)+" -gapext = "+str(ge)+" -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output = io.readFASTA("temp/tcoffee-out.fasta")
	return output

def createScoreMatrixFromAlignment(filename, output):
	""" Generate substitution matrix from a file containing alignment

	Args:
		filename (str): path to a fasta file containing aligned sequences
		output (str): path to a pickle file for the matrix to be stored

	Returns:
		my_lom (array): log odds matrix


	"""

	c_align = AlignIO.read(filename, "fasta")
	summary_align = AlignInfo.SummaryInfo(c_align)
	replace_info = summary_align.replacement_dictionary(["*"])
	my_arm = SubsMat.SeqMat(replace_info)
	
	# Add pseudocounts
	for m in my_arm:
		my_arm[m]+=1
	my_lom = SubsMat.make_log_odds_matrix(my_arm)
	
	pickle.dump( my_lom, open( output, "wb" ) )
	return my_lom

def tcoffeeMakeMatrix(M, output, tc_mismatch):
	""" Convert a numpy array corresponding to a substitution matrix to T-COFFEE format

	T-COFFEE needs a integer-based matrix in a special formatting

	Args:
		M (array): numpy substitution matrix
		output (str): path to file for saving the matrix
		tc_mismatch (int): default mismatch penalty 


	"""
	
	alpha = BIOALPHABET[:]
	string = "   "
	for i in range(0,len(alpha)-1):	
		string = string + alpha[i]+"  "
	string = string + alpha[i+1]+"  *\n"
	
	r = range(0,len(M))
	for i in r:
		string = string+alpha[r.index(i)]
		for j in r:
			val = str(int(round(max([M[i][j],-20])))) 
			if len(val)==1: string = string+"  "+val+""
			if len(val)==2: string = string+" "+val+""
			if len(val)==3: string = string+" "+val+""
		string = string+" "+str(tc_mismatch)+" \n"
	string = string+"*"
	
	for i in r:
		string = string+" "+str(tc_mismatch)+""
	string = string+" "+str(tc_mismatch)

	fo = open(output, "w")
	fo.write(string)
	fo.close()
