import string
import numpy as np
from Bio import SeqIO
import time
from types import *
import os

#redefined functions for handling FASTA
def exportGroupOfSequencesToFasta(num_seq_MSA,output):
	f = open(output, 'w')
	f.close
	f = open(output, 'a')
	for i,seq in enumerat(num_seq_MSA):
		f.write(">"+str(i)+"\n")
		if type(seq)==ListType:
			f.write("".join(seq))
		if type(seq)==StringType or type(seq)==UnicodeType:
			f.write(seq)
		f.write("\n")
	f.close

def readFASTA(file):
	newS=[]
	for seq_record in SeqIO.parse(file, "fasta"):
		newS.append(list(seq_record.seq))
	return newS



#A series of functions for performing multiple alignment using T-COFFEE
def tcoffeeCombineAlignments(sequences,alignedFASTAfiles):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-inx.fasta")
	string_of_fastas=""
	for file in alignedFASTAfiles:
		string_of_fastas+=file+" "
	command="./t_coffee 'temp/tcoffee-inx.fasta' -aln "+string_of_fastas+" -quiet -output fasta -outfile 'temp/tcoffee-outx.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-outx.fasta")
	return output


def tcoffeeAlignment(sequences,go,ge):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command="./t_coffee 'temp/tcoffee-in.fasta' -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeAlignmentWithID(sequences,go,ge):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command="./t_coffee 'temp/tcoffee-in.fasta' -matrix ID  -gapopen="+str(go)+" -gapext="+str(ge)+" -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeAlignmentWithTemplates(sequences,go,ge,matrixfile):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command="./t_coffee 'temp/tcoffee-in.fasta' -template_file 'template_file.tmpl'  -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeAlignmentWithMatrix(sequences,go,ge,matrixfile):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command="./t_coffee 'temp/tcoffee-in.fasta' -in=X"+matrixfile+"  -gapopen="+str(go)+" -gapext="+str(ge)+" -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-out.fasta")
	return output

def tcoffeeFileAlignment(file,go,ge):
	exportGroupOfSequencesToFasta(sequences,"temp/tcoffee-in.fasta")
	command="./t_coffee '"+file+"' -matrix ID  -gapopen="+str(go)+" -gapext="+str(ge)+" -quiet -output fasta -outfile 'temp/tcoffee-out.fasta'"
	os.system(command)
	output=readFASTA("temp/tcoffee-out.fasta")
	return output





