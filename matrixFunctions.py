from __future__ import division

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import operator
import os
import pickle
import random
import string

from Bio import AlignIO
from Bio import Alphabet
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC
from os import listdir
from os.path import isfile
from os.path import join
from types import *

BIOALPHABET = "ARNDCQEGHILKMFPSTWYVBZX".split() #Protein alphabet

#generate a substitution matrix from a fasta input file
def createScoreMatrixFromAlignment(filename,output,print_=False):
	c_align = AlignIO.read(filename, "fasta")
	summary_align = AlignInfo.SummaryInfo(c_align)
	replace_info = summary_align.replacement_dictionary(["*"])
	my_arm = SubsMat.SeqMat(replace_info)
	
	#add pseudocounts
	for m in my_arm:
		my_arm[m]+=1
	my_lom = SubsMat.make_log_odds_matrix(my_arm)
	
	pickle.dump( my_lom, open( output, "wb" ) )
	return my_lom


#convert a numpy array corresponding to a substitution matrix to T-COFFEE format
def tcoffeeMakeMatrix(M,alpha,output,tc_mismatch,print_=False):
	print "Creating the T-coffee matrix..."
	
	alpha=BIOALPHABET[:]
	string="   "
	for i in range(0,len(alpha)-1):	
		string=string+alpha[i]+"  "
	string=string+alpha[i+1]+"  *\n"
	
	r=range(0,len(M))
	if print_:
		print r
	
	for i in r:
		string=string+alpha[r.index(i)]
		for j in r:
			val=str(int(round(max([M[i][j],-20])))) 
			if len(val)==1: string=string+"  "+val+""
			if len(val)==2: string=string+" "+val+""
			if len(val)==3: string=string+" "+val+""
		string=string+" "+str(tc_mismatch)+" \n"
	string=string+"*"
	
	for i in r:
		string=string+" "+str(tc_mismatch)+""
	string=string+" "+str(tc_mismatch)

	fo = open(output, "w")
	fo.write(string)
	fo.close()
	if print_: #show on screen
		print string


