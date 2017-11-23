from __future__ import division
import os
from os import listdir
from os.path import isfile, isdir, join
import time
import re
import json
import math
import numpy as np
import string
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from types import *
import collections
from termcolor import colored

#Read and write json,FASTA,txt file functions
def readTxtFile(file):
	f = open(file, 'r')
	text=f.read()
	f.close()
	return text

def filesInPath(mypath,ext=[]):
	if len(ext)==0:
		onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f!=".DS_Store" ]
	else :
		onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f!=".DS_Store" and f[-len(ext):]==ext]
	return onlyfiles

def dirsInPath(mypath):
	return [ f for f in listdir(mypath) if isdir(join(mypath,f)) ]
	
def writeTxtFile(file,string):
	f = open(file, 'w')
	f.write(string)
	f.close()

def readJsonFile(file):
	with open(file) as data_file:
		data = json.load(data_file)
	return data

def readFASTA(file):
	return [list(seq_record.seq) for seq_record in SeqIO.parse(file, "fasta")]
		
def readFASTAandIDs(file):
	newS=[]
	ids=[]
	for seq_record in SeqIO.parse(file, "fasta"):
		ids.append(seq_record.id)
		newS.append(list(seq_record.seq))
	return newS,ids

def exportGroupOfSequencesToFASTA(num_seq_MSA,output):
	f = open(output, 'w')
	f.close
	f = open(output, 'a')
	for i,seq in enumerate(num_seq_MSA):
		f.write(">"+str(i)+"\n")
		if type(seq)==ListType:
			f.write("".join(seq))
		if type(seq)==StringType or type(seq)==UnicodeType:
			f.write(seq)
		f.write("\n")
	
	f.close

def exportGroupOfSequencesToFASTAwithIDs(num_seq_MSA,ids,output):
	f = open(output, 'w')
	f.close
	f = open(output, 'a')
	for i,seq in enumerate(num_seq_MSA):
		f.write(">"+ids[i]+"\n")
		if type(seq)==ListType:
			f.write("".join(seq))
		if type(seq)==StringType or type(seq)==UnicodeType:
			f.write(seq)
		f.write("\n")
	f.close

