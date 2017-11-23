from __future__ import division
import os
import sys, getopt
import time
import numpy as np
import random
from random import shuffle
import math
from types import *
import string

#global
hmmer_folder="HMMER/binaries-"
hmmer_profiles="HMMER/profiles"
hmmer_databases="HMMER/databases"
hmmer_options=" --music --wnone --enone  --plaplace " # --singlemx --mxfile 'id.matrix' "
hmmer_threshold=0.0
hmmer_correctLabel="!CORRECT!"

#jackhmmer search
def hmmJackSearch(query,id,targetFasta,toReturn,dataset="3"):
	fun.exportGroupOfSequencesToFASTAwithIDs([query],[id],"temp/jackhmmerinput.fasta")
	c_line=hmmer_folder+dataset+"/jackhmmer  -N 3 --fragthresh 1   --wnone --enone  --plaplace  --mxfile 'id.matrix' --max  --noali --domZ "+str(toReturn)+" -Z "+str(toReturn)+"  --tblout temp/jackhmmer.out  'temp/jackhmmerinput.fasta' "+targetFasta +" >temp/verbose.txt"
	os.system(c_line)
	rank=readHMMsearchOutput("temp/jackhmmer.out")
	return rank

def hmmSearch(id,targetFasta,toReturn,dataset="3"):
	profile=hmmer_profiles+"/"+id+".hmm"
	c_line=hmmer_folder+dataset+"/hmmsearch  --max -E 1e+20000000 --domE 1e+20000000  --noali --nonull2 --domZ "+str(toReturn)+" -Z "+str(toReturn)+" --tblout temp/hmmsearch.out  "+profile+" "+targetFasta+" >temp/verbose.txt"
	os.system(c_line)
	rank=readHMMsearchOutput("temp/hmmsearch.out")
	return rank

def hmmScan(sequence,id,toReturn=1000,dataset="3"):
	fun.exportGroupOfSequencesToFASTAwithIDs([sequence],[id],"temp/hmmscaninput.fasta")
	c_line=hmmer_folder+dataset+"/hmmscan --max -E 1e+20000000 --domE 1e+20000000  --nonull2 --noali --domZ "+str(toReturn)+" --tblout temp/hmmscan.out  '"+hmmer_databases+"/database' temp/hmmscaninput.fasta"+" >temp/verbose.txt"
	os.system(c_line)
	rank=readHMMsearchOutput("temp/hmmscan.out")
	return rank

#read and convert hmmer outpur
def readHMMsearchOutput(file):
	f = open(file, 'r')
	top=[] #in case noone passes the threshold
	scores=[]
	for line in f:
		if line[0]!="#":
			#print line.split()
			#line.split()[0].split("-")[1].split(".")[0]
			# print top
			# time.sleep(1000)
			top.append(line.split()[0])
			scores.append(float(line.split()[5]))
	f.close()
	#print len(top)
	return top,scores

#build hmmer profile from folder
def createHMMprofiles(folder,threshold,subfolder,dataset):
	cliques=fun.filesInPath(folder)
	for clique in cliques:
		input_file='"'+folder+'/'+clique+'"'
		c_line=hmmer_folder+dataset+'/hmmbuild --fragthresh 1 --symfrac '+str(threshold)+hmmer_options+' "'+hmmer_profiles+'/'+subfolder+'/'+clique+'.hmm" '+input_file
		os.system(c_line)
	return 1

#build profile database
def createHMMDatabase(dataset="3"):
	#create
	print "Create database from HMMs..."
	print "Removing old..."
	os.system("rm -Rf "+hmmer_databases+"/database")
	os.system("rm -Rf "+hmmer_databases+"/database.h3f")
	os.system("rm -Rf "+hmmer_databases+"/database.h3i")
	os.system("rm -Rf "+hmmer_databases+"/database.h3m")
	os.system("rm -Rf "+hmmer_databases+"/database.h3p")
	print "Concatenating..."
	c_line="cat "
	cliques=fun.filesInPath(hmmer_profiles)
	for clique in cliques:
		c_line=c_line+' "'+hmmer_profiles+"/"+clique+'" '
	c_line=c_line+' > "'+hmmer_databases+'/database"'
	
	#run the concat command
	os.system(c_line)
	
	#compress
	print "Compressing database..."
	c_line=hmmer_folder+dataset+"/hmmpress "+hmmer_databases+"/database"
	os.system(c_line)
	return 1

#build hmmer profile from alignment
def createHMMprofileFromAlignment(alignment,name,dataset="3"):
	os.system("rm -Rf "+hmmer_profiles+"/"+name+'.hmm') #delete previous correct hmm
	fun.exportGroupOfSequencesToFasta(alignment,"temp/"+name+".fasta")
	input_file="temp/"+name+".fasta"
	c_line=hmmer_folder+dataset+'/hmmbuild -o hmm.out  --symfrac '+str(hmmer_threshold)+hmmer_options+' "'+hmmer_profiles+'/'+name+'.hmm" '+input_file
	os.system(c_line)
	os.system("rm -Rf 'temp/"+name+".fasta'")









