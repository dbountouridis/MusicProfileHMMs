## -*- coding: cp1252 -*-
from __future__ import division
import os
from os import listdir
from os.path import isfile, isdir, join
import time
import re
import json
import pylev
import math
import numpy as np
from random import shuffle
import pickle
from scipy import stats
from scipy import spatial
import math
import string
import math
import random
from types import *
import functions as fun
import networkx as nx
import string
import matplotlib.pyplot as plt
import plots as myplots
import tcoffee as tcof
from sklearn.metrics import f1_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
#import smirnov_grubbs as grubbs
#import outlier as out
import seaborn as sns
import pandas as pd
import matrix as matrx
import hmmProf
import collections
import entropy
import networkx
from ballotbox.ballot import BallotBox
from ballotbox.singlewinner.preferential import BordaVoting
import rank
import rank_aggregators as r
from sets import Set
import shutil
from itertools import permutations
import hmmer


a="A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X "
Alphabet=a.split()

def f(sizes,scores,labels,ids):
	sizes=np.asarray(sizes)
	scores=np.asarray(scores)
	labels=np.asarray(labels)
	ids=np.asarray(ids)
	index=np.argsort(scores)
	scores=scores[index]
	labels=labels[index]
	sizes=sizes[index]
	ids=ids[index]
	return sizes,scores,labels,ids

def removeGapColumns(sequences):
	for i,seq in enumerate(sequences):
		if i==0:U=Set(np.where(np.array(list(seq))=="-")[0])
		else:
			s=Set(np.where(np.array(list(seq))=="-")[0]) 
			U=U&s
	U=list(U)
	U.sort()
	U.reverse()
	ns=[]
	for sequ in sequences:
		seq=list(sequ)
		#print len(seq)
		for index in U:	seq.pop(index)
		ns.append("".join(seq))
	return ns

def makefolder(folder):
	if not os.path.exists(folder): os.makedirs(folder)

def makeAlignmentPrototypes(aligned,type="fusion",k=5):
	if type=="fusion":
		fu=fun.fusionConsensus(aligned,"temp/fin","temp/fout")
		fu=fu.replace(",","")
		return fu

	if type =="mj":
		mj=fun.consensusFromSequencesCnf(aligned,k/10,"-",0)
		#print mj
		return str(mj)
	
	if type=="random":
		rc=fun.randomConsensus(aligned)
		return "".join(rc)

def sequenceToMSAviaSplitting(seq,split="half",depth=0.2,firstSplitAt=-1):
	S=[]
	divi=[]
	S.append(seq)
	l=len(seq)
	divi.append(0)
	firstSplit=False

	if split in ["half","random","fuse"]:
		for j,q in enumerate(S):
			#if the sequence has not been split before and its length is larger than t% of the original
			if divi[j]==0 and len(q)>=4 :
				if split in ["half","fuse"]:
					s1=q[0:int(len(q)/2)]
					s2=q[int(len(q)/2):]

				#randomly split the sequence, 70-30% is the maximum split
				if split=="random":
					rand=min(max(0.3,random.random()),0.7)
					s1=q[0:int(len(q)*rand)]
					s2=q[int(len(q)*(1-rand)):]

				if firstSplitAt>0 and firstSplit==False:
					s1=q[0:firstSplitAt]
					s2=q[firstSplitAt:]
					firstSplit=True
				
				#if the resulting splits have absolute length larger than
				if len(s1)>=2 and len(s2)>=2:
					divi[j]==1
					S.append(s1)
					S.append(s2)
					divi.append(0)
					divi.append(0)
	if split in ["frans","fuse"]:
		S=S[:3]
		for k in range(1,4):
			s=np.roll(np.array(seq), int(k/4*len(seq)))
			S.append(s.tolist())

	if split=="empty":
		for k in range(0,4): S.append(seq)

	if split=="permutations":
		S=[]
		lenAlphabet=3
		alpha="ABCDEFGHIJKLMNOPQRST"
		alpha=list(alpha)
		dict_={}
		#print seq,len(seq)
		for i in range(lenAlphabet):			
			part="".join(seq[int(len(seq)/lenAlphabet*i):int(len(seq)/lenAlphabet*(i+1))])
			dict_.update({alpha[i]: "".join(part) })
		

		keywords = [''.join(i) for i in permutations(alpha[:lenAlphabet],  lenAlphabet)]
		for k in keywords:
			#print k
			s=""
			for char in list(k): s+=dict_[char]
			S.append(list(s))
	

	
	if split=="four":
		for k in xrange(2):
			rand=min(max(0.3,random.random()),0.7)
			s1=seq[0:int(len(seq)*rand)]
			s2=seq[int(len(seq)*(1-rand)):]
			S.append(s1)
			S.append(s2)

	# while len(S)>4:
	# 	S.pop(-1)



	F=[]
	for j,q in enumerate(S):
		if (  len(F)%2==1 and len(F)>=9): continue #j/len(S)>depth
		#print j,"".join(q)
		s=[]
		while len(s)<len(seq):
			s+=q
		#s=s[:len(seq)]
		F.append(s)
	return F

def computeNullEmissionProbs(alignment):
	#create null Emission probs for insert states from the data
	EmissionProbs={}
	for char in Alphabet:
		EmissionProbs.update({char:0})
	string=""
	#Emission counts for hmm

	for s in alignment:
		string+="".join(s).replace("-","")
	counts=collections.Counter(string)
	for char in counts:
		if char!="-": EmissionProbs[char]+=counts[char]
	sum=0
	for char in Alphabet: sum+=EmissionProbs[char]
	for char in Alphabet: EmissionProbs[char]= EmissionProbs[char]/sum
	print "nullEmission probs:",EmissionProbs
	
	todel=[]
	for char in Alphabet: 
		if EmissionProbs[char]==0: todel.append(char)
	EmissionProbsBefore=EmissionProbs.copy()
	for char in todel: 
		del EmissionProbs[char]
	newAlphabet=[]
	for char in Alphabet:
		if char not in todel: newAlphabet.append(char)

	print "nullEmission probs (deleted empty):", EmissionProbs
	print "Alphabet:",newAlphabet
	
	return EmissionProbs,newAlphabet,EmissionProbsBefore

def createHMMmodel(aligned,EmissionProbs,Alphabet,d_inertia=0.5, e_inertia=0.5,s_th=0.001,sw_=1,pct_=1,pec_=1,grms_=0.39):
	
	
	print fun.printMSA(aligned)
	#currentEmissionProbs,currentAlphabet,_=computeNullEmissionProbs(aligned)

	if len(aligned)<=4:
		grms=0.99
	else:
		grms=grms_#0.5+random.random()*0.49

	pct=pct_#random.random()*10
	sw=sw_#int(10+random.random()*400)
	pec=pec_#random.random()*10

	model,blueprint=hmmProf.profileHMM(aligned,Alphabet,gapRatiotoMatchStates=grms,pseudoCountTransition=pct,sequenceWeight=sw,pseudoEmissionCount=pec,plot_=False,uniformInsertionProbs=False,nullEmission=EmissionProbs)
	str1 = ''.join(str(e) for e in blueprint.tolist())
	print str1, "<- Insert states"

	

	T=[]
	for t in aligned: # we train only on the original sequence! IMPORTANT
		t="".join(t)
		for b in range(0,1):T.append(list(t.replace("-","")))
	model.train(T,max_iterations=500,distribution_inertia=d_inertia, edge_inertia=e_inertia,stop_threshold=s_th,algorithm='baum-welocalh')

	return model
		
def analyseMSA(aligned):
	numOfsequences=len(aligned)
	lengthOfMSA=len(aligned[0])
	a=np.array(aligned)
	fun.printMSA(aligned)

	d=[]
	ratios=[]
	percentageOfgaps=[]
	for i in xrange(10):
		d.append(list(" "*lengthOfMSA))
	#print d

	for i in xrange(lengthOfMSA):
		column=a[:,i]
		#print column[np.where(column!="-")]
		counts=collections.Counter(column[np.where(column!="-")])
		percentageOfgaps.append(len(np.where(column=="-")[0])/numOfsequences)
		mx=0
		for char in counts:
			freq=counts[char]
			if freq>mx: mx=freq	
		ratio=mx/numOfsequences
		#print ratio
		ratios.append(ratio)
		for j in xrange(int(math.floor(ratio*10))):
			d[j][i]="|"
	for i in range(len(d)-1,0,-1):
		print "".join(d[i])
	
	#print numOfsequences,lengthOfMSA,np.mean(ratio),np.mean(percentageOfgaps)
	return numOfsequences,lengthOfMSA,np.mean(ratios),np.mean(percentageOfgaps),ratios

def randomSustrings(l, k, m, result=[]):
    if len(result) == m or len(l) == 0:
        return result
    else:
        if isinstance(l, str):
            l = [l]
        part_num = random.randint(0, len(l)-1)
        partition = l[part_num]
        start = random.randint(0, len(partition)-k)
        result.append(partition[start:start+k])
        l.remove(partition)
        l.extend([partition[:start], partition[start+k:]])
        return randomSustrings([part for part in l if len(part) >= k], k, m, result)

def shuffleBoth(list1,list2):
	l1=list1[:]
	shuffle(l1)
	l2=[]
	for id in l1: l2.append(list2[list1.index(id)])
	return l1,l2

def degradeSequence(q_sequence,degradationType="half",folderForAdditions=[],queryIndex=0):
	l=len(q_sequence)
	if degradationType=="normal": return q_sequence
	if degradationType=="half":
		index=int(random.random()*(l/2))
		new=q_sequence[index:index+int(l/2)]
		return new
	if degradationType=="double":
		new=q_sequence[:]
		new+=q_sequence
		return new
	if degradationType=="additions":
		files=fun.filesInPath(folderForAdditions)
		print "len:",len(files)
		#select a random sequence
		pick=int(random.random()*len(files))
		while pick==queryIndex:
			pick=int(random.random()*len(files))
			print pick
		sequencesA,idsA=fun.readFASTAandIDs(folderForAdditions+"/"+files[pick])
		shuffle(sequencesA)
		addedSequence=sequencesA[0]
		index=int(random.random()*min(l*0.66,len(addedSequence)))
		a_sequence=addedSequence[index:min(index+int(l*0.33),index+len(addedSequence))]
		q_sequence="".join(q_sequence).replace("-","")
		a_sequence="".join(addedSequence).replace("-","")
		#where to place
		index=int(random.random()*l)
		new=q_sequence[:index]+a_sequence+q_sequence[index:]
		return new
	if degradationType=="randomChars":
		#probable chars
		alpha=list(set(q_sequence))
		new=q_sequence[:]
		for k in range(0,l):
			if random.random()<0.25: new[k]=alpha[int(random.random()*len(alpha))]
		return new
	if degradationType=="permutations":
		r=randomSustrings("".join(q_sequence), int(len(q_sequence)*0.2), 5,result=[])
		shuffle(r)
		new=""
		for rr in r: new+=rr
		return new
	if degradationType=="repetitions":
		newstring=[]
		window=4
		for index in range(0,l,window):
			s=q_sequence[index:index+window]
			if random.random()>0.5: #normal case
				newstring+=s
			else:
				if random.random()>0.5: #repeat x2
					newstring+=s
					newstring+=s
		return newstring
def expandSplits(splits,L):
	newSplits=[]
	for s in splits:
		while len(s)<L:
			s+=s
		newSplits.append(s[:len(splits[0])])
	return newSplits
# halfing the sequences until the have subsequences are of length 4
def halfSplit(seq):
	seeds=[]
	seeds.append(seq)
	splits=[]
	splits.append(seq)

	while len(seeds)>=1:
		#print seeds
		#time.sleep(1)
		s1=seeds[0][0:int(len(seeds[0])/2)]
		s2=seeds[0][int(len(seeds[0])/2):]
		if len(s1)>=4 and len(s2)>=4:
			seeds.append(s1)
			seeds.append(s2)
			splits.append(s1)
			splits.append(s2)
		seeds.pop(0)
	s2=[]
	for s in splits:
		s2.append(s[:len(seq)])
	splits=s2[:]
	return splits
def shiftSplit(seq,maximumNumber):
	seeds=[]
	seeds.append(seq)
	for k in range(1,maximumNumber):
		s=np.roll(np.array(seq), int(k/maximumNumber*len(seq)))
		seeds.append(s.tolist())
	return seeds
def permutationsSplit(seq):
	lenAlphabet=4  # this means there are 24 permutations
	alpha="ABCDEFGHIJKLMNOPQRST"
	alpha=list(alpha)
	dict_={}
	seeds=[]
	#Create dictionary
	for i in range(lenAlphabet):			
		part=seq[int(len(seq)/lenAlphabet*i):int(len(seq)/lenAlphabet*(i+1))]
		dict_.update({alpha[i]: part })
	keywords = [''.join(i) for i in permutations(alpha[:lenAlphabet],  lenAlphabet)]
	for k in keywords:
		s=[]
		for char in list(k): s+=dict_[char]
		seeds.append(list(s))
	return seeds		
def plotMatrix(matrix, filename):
	alpha=[]
	for key in matrix.keys():
		alpha.append(key[0])
		alpha.append(key[1])
	alpha=np.sort(list(set(alpha))).tolist()

	Alphabet = alpha
	M = np.zeros((len(Alphabet), len(Alphabet)))
	for i, char1 in enumerate(Alphabet):
		for j, char2 in enumerate(Alphabet):
			if (char1, char2) in matrix.keys():
				M[i][j] = matrix[(char1, char2)]
			else:
				M[i][j] = matrix[(char2, char1)]
	corr = pd.DataFrame(M, index=alpha, columns=alpha)
	print corr

	sns.set_context("notebook", font_scale=0.8, rc={"lines.linewidth": 1.0})
	sns.set_style({'font.family': 'serif', 'font.serif': ['Times New Roman']})
	#sns.set_style({'text.usetex': True})

	# Generate a mask for the upper triangle
	mask = np.zeros_like(corr, dtype=np.bool)
	mask[np.triu_indices_from(mask)] = True
	# Set up the matplotlib figure
	f, ax = plt.subplots(figsize=(11, 9))

	# Generate a custom diverging colormap
	cmap = sns.diverging_palette(220, 10, as_cmap=True)

	# Draw the heatmap with the mask and correct aspect ratio
	sns.heatmap(corr, cmap=cmap,  cbar_kws={"shrink": 1.0}, ax=ax, annot=True,fmt=".1f")
	#sns.plt.savefig("plots/" + filename + ".pdf", format="pdf")
	
	#plt.show()
def makeArtificialVersionsFromSequence(seq,split="half",maximumNumber=6):
	S=[]
	divI=[]
	S.append(seq)
	indeces=range(len(seq))
	I=[]
	I.append(indeces)
	if split=="nothing":
		splits=[]
		indexSplits=[]
		for i in range(maximumNumber):
			splits.append(seq)
			indexSplits.append(indeces)
		return splits,indexSplits
	if split=="random":
		splits=[]
		indexSplits=[]
		for i in range(maximumNumber):
			r=range(len(seq))
			shuffle(r)
			splits.append(np.array(seq)[r].tolist())
			indexSplits.append(np.array(indeces)[r].tolist())
		return splits,indexSplits
	if split=="half":
		splits=expandSplits(halfSplit(seq),len(seq))
		indexSplits=expandSplits(halfSplit(indeces),len(indeces))
		return splits[:maximumNumber],indexSplits[:maximumNumber]
	if split=="shift":
		splits=shiftSplit(seq,maximumNumber)
		indexSplits=shiftSplit(indeces,maximumNumber)
		#print splits,indexSplits
		return splits[:maximumNumber],indexSplits[:maximumNumber]
	if split=="halvesAndShift":
		#halves
		splits=expandSplits(halfSplit(seq),len(seq))
		indexSplits=expandSplits(halfSplit(indeces),len(indeces))
		#shifts
		splits2=shiftSplit(seq,maximumNumber)
		indexSplits2=shiftSplit(indeces,maximumNumber)
		#combine 3 and 3
		splits=splits[:int(maximumNumber/2)]+splits2[int(maximumNumber/2):maximumNumber]
		indexSplits=indexSplits[:int(maximumNumber/2)]+indexSplits2[int(maximumNumber/2):maximumNumber]
		return splits[:maximumNumber],indexSplits[:maximumNumber]
	if split=="permutations":
		splits= permutationsSplit(seq)
		indexSplits= permutationsSplit(indeces)
		return splits[:maximumNumber],indexSplits[:maximumNumber]

###################### Instructions ##############################
# first run the createdataset function for CSVproperSmall and NLBproperSmall(smaller alphabet)
# then run the matlab function experiment.m
# then the experiment function here
# then the analysis



a="A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X "
Alphabet=a.split()

dataset="NLBproperSmall"



create_datasets=0
experiment=1
frequenciesForHmmer=0
multiplication=40
analysis=0


#type is the type of degradataion applied to the query
#type2 is the degradation applied to the sequences in the family


# altprot_file="temp/altProt-"+dataset+".pkl"
# dataset_file="temp/dataset-"+dataset+".pkl"
# alignments_file="temp/alignments-"+dataset+".pkl"
# output_aggr="dataframe_aggr_"+dataset+".pkl"
# output_aggr2="dataframe_aggr2_"+dataset+".pkl"
# hmm_file="temp/hmm-"+dataset+".pkl"
# hmmer_file="temp/hmmer-"+dataset+".pkl"
# jackhmmer_file="temp/jackhmmer-"+dataset+".pkl"
# hmmer_file_reverse="temp/hmmer-reverse-"+dataset+".pkl"
# hmm_file_inverse="temp/hmm-inverse-"+dataset+".pkl"
# msapath="MultipleAlignments/"+dataset+"/"
# msapathqueries="MultipleAlignments/"+dataset+"-targets/"
# mypath="NotAligned"
# alignments="MultipleAlignments/"
# dirs=fun.dirsInPath(mypath)
# hmms_folder="hmmFolder/"


#My profile HMM settings
di=0.5
ei=0.5
#pseudo=1#int(random.random()*10)+1
sw=10#int(random.random()*10)+1
pct=1
pec=1/12.
s_th=0.0000000001# pow(10,-int(random.random()*5+1))
grms=0.95#0.1+random.random()*0.8 #0.39


#Offline datasets creation for faster processing
def create_dataset(dataset="NLBproperSmall"):
	database="Notaligned/"+dataset
	files=fun.filesInPath(database)

	#create folder
	if not os.path.isdir("MultipleAlignments"): os.makedirs("MultipleAlignments")
	if not os.path.isdir("MultipleAlignments/prototypes"): os.makedirs("MultipleAlignments/prototypes")
	if not os.path.isdir("MultipleAlignments/prototypes/"+dataset): os.makedirs("MultipleAlignments/prototypes/"+dataset)
	if not os.path.isdir("MultipleAlignments/prototypes/"+dataset): os.makedirs("MultipleAlignments/queries/"+dataset)
	if not os.path.isdir("MultipleAlignments/queries/"+dataset): os.makedirs("MultipleAlignments/queries/"+dataset)
	for degradationType in ["additions","normal","repetitions","permutations","randomChars","double" ,"half"]:
		if not os.path.isdir("MultipleAlignments/prototypes/"+dataset+"/"+degradationType): 
			os.makedirs("MultipleAlignments/prototypes/"+dataset+"/"+degradationType)

	SEQ=[]
	IDS=[]
	#Create sequences , prototypes, degradations
	for j in range(0,len(files[:])):
		sequences,ids=fun.readFASTAandIDs("Notaligned/"+dataset+"/"+files[j])
		SEQ+=sequences
		IDS+=[files[j].replace(".fasta","")+"-"+str(i) for i in range(len(ids))]

	#degrade
	for degradationType in ["additions","normal","repetitions","permutations","randomChars","double" ,"half"]:
		SEQD=[]
		#create aligned prototypes
		for i in range(len(SEQ)):
			sequence=SEQ[i]
			id=IDS[i]
			print id,degradationType, i/len(SEQ)
			#degrade
			degraded=degradeSequence(sequence,degradationType=degradationType,folderForAdditions="Notaligned/"+dataset,queryIndex=i)
			print "".join(sequence)
			print degraded

			SEQD.append(degraded)

			#create versions
			versions,indeces=makeArtificialVersionsFromSequence(degraded,split="half",maximumNumber=6)

			fun.printMSA(versions)
			aligned=fun.runMafftAlignmentWithSettings(versions,2,1,method="localpair",allowshift=True)
			fun.printMSA(aligned)
			fun.exportGroupOfSequencesToFASTA(aligned,"MultipleAlignments/prototypes/"+dataset+"/"+degradationType+"/"+id)
		fun.exportGroupOfSequencesToFASTAwithIDs(SEQD,IDS,"MultipleAlignments/queries/"+dataset+"/"+degradationType)
						
	


if create_datasets:
	create_dataset(dataset=dataset)

if frequenciesForHmmer:
	a="A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X "
	database="Notaligned/"+dataset
	databaseFiles=fun.filesInPath(database)
	databaseCliques=[]
	for file in databaseFiles:
		sequences,_=fun.readFASTAandIDs(database+"/"+file)
		databaseCliques.append(sequences)
	#EMission probs from the whole set
	AllSeqs=[]
	for clique in databaseCliques: 
		for seq in clique: AllSeqs.append(seq)
	EmissionProbs,newAlphabet,_=computeNullEmissionProbs(AllSeqs)
	print EmissionProbs
	s=""
	for i,char in enumerate(a.split()):
		val=0
		if char in EmissionProbs.keys():
			val=EmissionProbs[char]
		s+="f[%d]=%.5f;\n"%(i,val)
	print s

def expriment(dataset="NLBproperSmall"):
	
	for queryDegradationType in ["normal"]:
		for degradationType in ["additions","normal","repetitions","permutations","randomChars","double" ,"half"]:

			#read the queries
			QuerySEQ,QueryIDS=fun.readFASTAandIDs("MultipleAlignments/queries/"+dataset+"/"+queryDegradationType)

			#read the target MSA database
			database="MultipleAlignments/prototypes/"+dataset+"/"+degradationType
			databaseFiles=fun.filesInPath(database)

			#read the target queries for pairwise alignment
			targetSEQ,targetIDS=fun.readFASTAandIDs("MultipleAlignments/queries/"+dataset+"/"+degradationType)

			pairwise=0
			if pairwise:
				print "Pairwise ", queryDegradationType, degradationType
				M=np.zeros((len(QuerySEQ),len(targetSEQ)))
				L=np.zeros((len(QuerySEQ),len(targetSEQ)))

				for i,query in enumerate(QuerySEQ):
					print i/len(QuerySEQ)
					query_id=QueryIDS[i].split("-")[0]
					for j,target in enumerate(targetSEQ):
						target_id=targetIDS[j].split("-")[0]

						score=fun.pairwiseBetweenSequences("".join(query),"".join(target), gapopen=-0.8,gapext=-0.5,type_="global")
						score=score/np.max([len(query),len(target)])
						M[i][j]=score
						L[i][j]=target_id==query_id

					pickle.dump((M,L),open("results/"+dataset+"/pairwise-"+queryDegradationType+"-"+degradationType+".pickle","wb" ) )
	



			hmmer_offline=1				
			if hmmer_offline:
				print "HMMER ", queryDegradationType, degradationType
				M=np.zeros((len(QuerySEQ),len(targetSEQ)))-1000000000000
				L=np.zeros((len(QuerySEQ),len(targetSEQ)))

				shutil.rmtree("HMMER/profiles")
				os.makedirs("HMMER/profiles")
				

				#read
				print "Create profiles..."
				databaseCliques=[]
				databaseCliqueIds=[]
				for file in databaseFiles:
					sequences,_=fun.readFASTAandIDs(database+"/"+file)
					databaseCliques.append(sequences)
					databaseCliqueIds.append(file)

					a=[]
					for k in range(0,multiplication):
						for seq in sequences: a.append(seq)
					aligned=a[:]
					hmmer.createHMMprofileFromAlignment(aligned,file,dataset)
				hmmer.createHMMDatabase(dataset)

				
				
				AP=[]
				for j,query in enumerate(QuerySEQ):
					query_id=QueryIDS[j]
					query_clique_id=QueryIDS[j].split("-")[0]
					print j/len(QuerySEQ)

					rank,scores=hmmer.hmmScan(query,query_clique_id,dataset=dataset)
						
					# #save to a file
					# for k,target_id in enumerate(rank): 
					# 	D_hmmer.append([ query_id, "hmmer",scores[k],query_clique_id==target_id,type,type2])
					
					#some realtime feedback for me
					if len(rank)>0:
						labels=[]
						for k,item in enumerate(rank):
							index=targetIDS.index(item)
							L[j][index]=query_clique_id==item.split("-")[0]
							M[j][index]=scores[k]

				pickle.dump((M,L),open("results/"+dataset+"/hmmer-"+queryDegradationType+"-"+degradationType+".pickle","wb" ) )
					


	
	# if hmmer_offline:
	# 	Data=pd.DataFrame(data=D_hmmer,columns=columns3)
	# 	Data.to_pickle(hmmer_out)
				


if experiment:
	expriment(dataset=dataset)

analysis2=1
if analysis2:
	dataset="NLBproperSmall"
	files=fun.filesInPath("results/"+dataset)
	for file in files:
		f=file.replace(".pickle","").split("-")
		method=f[0]
		queryDegradationType=f[1]
		targetDegradationType=f[2]
		print method,queryDegradationType,targetDegradationType
		(M,L)=pickle.load( open("results/"+dataset+"/"+file, "rb" ))

		AP=[]
		for j in range(M.shape[0]):
			scores=M[j][:]
			labels=L[j][:]
			AP.append(average_precision_score(labels,scores))
		print np.mean(AP)



def smallAnalysis(df,type1="",types2="",reverse=False,output=""):
	columns=["Degradation Type","Average Precision","Classification Accuracy","Protype/Profile Method"]
	methodss=["mj-0.3","mj-0.0","mj-0.6","fusion","hmm-matlab-reverse","hmmer","profile-HMM","random","PW-exhaustive-"]
	mapm_=["MjV-0.3","MjV-0.0","MjV-0.6","Data Fusion","profile-HMM (Matlab)","HMMER","profile-HMM (Krogh)","Random Sequence","Exhaustive pairwise"]
	hue_order=["Exhaustive pairwise","Random Sequence","MjV-0.3","MjV-0.0","MjV-0.6","Data Fusion","profile-HMM (Matlab)","profile-HMM (Krogh)","HMMER"]

	for file in files:
		if os.path.isfile(file):
			print file
			df_=pd.read_pickle(file)
			#print df_
			df=df.append(df_)

	map_=["Repetitions","Permutations","Random Noise","Infusions","Double size" ,"Half size"]

	D=[]
	for type2 in types2:
		if reverse:
			AlignmentData=df.loc[(df['type'] == type1) & (df['type2'] == type2) ]
		else:
			AlignmentData=df.loc[(df['type'] == type2) & (df['type2'] == type1) ]


		#print
		#print "General MAP results, per method "+MSAmethod+" : ("+type1+","+type2+"):"
		
		# Hmmer doesnt return a full list of rankings. So we deal with that
		maps={}
		clas={}
		mins={}
		for query, group in AlignmentData.groupby("query"):
			ma_=100000
			for method, g in group.groupby("method"):
				if len(g)<ma_: 
					ma_=len(g)
			mins.update({query:ma_})

		g=[]
		for key in mins.keys():
			g.append(mins[key])
		#print "Average length of ranked list:",np.mean(g)
		
		query_name=[]
		method_name=[]
		for method, group in AlignmentData.groupby("method"):
			maps.update({method:[]})
			clas.update({method:[]})
			query_name=[]
			method_name.append(method)
			#print method
			for query, g in group.groupby("query"):
				if mins[query]>1:
					index=np.argsort(g["score"])[::-1]
					labels=np.array(g["class"])[index[0:mins[query]]].astype(int)
					scores=np.array(g["score"])[index[0:mins[query]]].astype(float)
					#print len(labels),len(scores)
					ap=average_precision_score(labels,scores)
					if math.isnan(ap): ap=0
					maps[method].append(ap)
					clas[method].append(int(labels[0]==1))
					query_name.append(query)
		
		#text=""
		for method in maps:
			#print method,type2,"  mean (std):",np.mean(maps[method]),"(",np.std(maps[method]),")"
			#print method,type2,"c mean (std):",np.mean(clas[method]),"(",np.std(clas[method]),")"
			for j,ap in enumerate(maps[method]):
				D.append([map_[types2.index(type2)],ap,clas[method][j], mapm_[methodss.index(method)]])
				D.append(["Average",ap,clas[method][j],mapm_[methodss.index(method)]])
		#print text

	

	#final figure
	Data=pd.DataFrame(data=D,columns=columns)

	#latex
	s=" & "
	grouped=Data.groupby("Degradation Type")
	for dtype, g in grouped:
		s+=dtype+" & "
	s+="\\\\ \n"
	grouped=Data.groupby("Protype/Profile Method")
	for method,group in grouped:
		s+=method+" & "
		grouped2=group.groupby("Degradation Type")
		for dtype, g in grouped2:
			ap=np.mean(g["Average Precision"])
			cl=np.mean(g["Classification Accuracy"])
			s+="%.2f & %.2f & "%(ap,cl)
		s+=" \\\\ \n"
	print s
	
	
	
	
	
	for evaluation in ["Average Precision","Classification Accuracy"]:
		#fig, ax = plt.subplots()
		# the size of A4 paper
		#fig.set_size_inches(20, 16)
		sns.set_style('whitegrid')
		sns.set_context("paper", font_scale=3)
		sns.set_style({'font.family':'serif', 'font.serif':['Times New Roman']})

		#	hue_order=["Exhaustive pairwise","Random Sequence","MjV-0.0","MjV-0.3","MjV-0.6","Data Fusion","profile-HMM (Matlab)","profile-HMM (Krogh)","HMMER"]
		colors = ["coral" ,"amber", "light blue","sky blue","dark sky blue","rose", "sage","grey green","moss green",  "blue grey","grey","lavender","lavender","lavender","light grey","amber","amber", "amber","amber","sage green", "green blue","amber","sage green", "green blue","amber","sage green", "green blue","amber"]
		sns.set_palette(sns.xkcd_palette(colors))
		f, ax = plt.subplots(figsize=(25, 6))
		g=sns.barplot(x="Degradation Type", y=evaluation, hue="Protype/Profile Method", data=Data,capsize=.0,order=["Double size" ,"Half size","Repetitions","Permutations","Random Noise","Infusions","Average"],hue_order=hue_order,)
		#sns.despine(left=True)
		#plt.legend(loc='upper right')
		box = ax.get_position() # get position of figure
		ax.set_position([box.x0, box.y0, box.width * 0.85, box.height]) # resize position
		ax.set(ylabel=evaluation, xlabel='')
		plt.legend(bbox_to_anchor=(1.00, 1), loc=2, borderaxespad=0.)
		#ax.set_title(MSAmethod)
		#sns.despine(offset=100,top=0, right=0, left=0, bottom=0)
		plt.setp(plt.gca().get_legend().get_texts(), fontsize=str(18))
		#plt.show()
		plt.savefig("plots/prototypeRobustness-"+output+"-"+evaluation+".pdf",format="pdf")
		plt.close()

if analysis:
	

	for MSAmethod in ["globalpairFalse"]:
		files=["ResultsProt/"+dataset+"-"+MSAmethod+"-hmmer","ResultsProt/"+dataset+"-"+MSAmethod+"-myhmm","ResultsProt/"+dataset+"-"+MSAmethod+"-pw","ResultsProt/"+dataset+"-"+MSAmethod+"-alt","ResultsProt/"+dataset+"-"+MSAmethod+"-mat"]
		columns3=["query","method","score","class","type","type2"]
		df=pd.DataFrame(columns=columns3)


		type1= "normal"
		types2=["normal"]
		smallAnalysis(df,type1,types2,output=dataset+"-normal2")
		smallAnalysis(df,type1,types2,reverse=True,output=dataset+"-reverse2")

		type1= "normal"
		types2=["repetitions","permutations","randomChars","additions","double" ,"half"]
		smallAnalysis(df,type1,types2,output=dataset+"-normal")
		smallAnalysis(df,type1,types2,reverse=True,output=dataset+"-reverse")

		









