""" A collection of reading/writing file functions. 

The collection includes reading (and/or) writing MIDI
FASTA, json and other file formats. Most of them are 
self-explenatory.

Example:
	See example.py


"""

from __future__ import division

import collections
import json
import math
import midi
import numpy as np
import os
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
from os import listdir
from os.path import isdir
from os.path import isfile
from os.path import join
from termcolor import colored
from types import *

__author__ = 'Dimitrios  Bountouridis'

def readMIDIfile(file):
	""" Read midi file.

	This function tries two different midi reading functions (midi type 1 and 2).

	Args:
		file (str): a file path

	Returns:
		Notes (dict): the pitch of notes for each chanel
		Durations (dict): the durations of notes for each chanel
		Onsets (dict): the notes' onset times
		minduration (float): the minimum note duration
		PitchIntervals (dict): the pitch interval of notes for each chanel

	
	"""


	try:
		print "Trying MIDI type 1..."
		Notes, Durations, Onsets, minduration, _ = readMIDItype1(file)
	except RuntimeError:
		print "Trying a different MIDI file format..."
		Notes, Durations, Onsets, minduration,_ = readMIDItype2(file)
	
	# Pitch intervals is a an essential sequence representation
	PitchIntervals = {}
	for channel in Notes.keys(): PitchIntervals.update({channel:[Notes[channel][i+1]-Notes[channel][i] for i in range(len(Notes[channel])-2)]})
	
	return Notes, Durations, Onsets, minduration, PitchIntervals

def readMIDItype1(file):
	pattern = midi.read_midifile(file)
	resolution = pattern._EventStream__resolution
	events = pattern.trackpool
	# We assume a 4 channel MIDI file
	Notes = {0:[],1:[],2:[],3:[]}
	Durations = {0:[],1:[],2:[],3:[]}
	Onsets = {0:[],1:[],2:[],3:[]}
	minduration = 10000000
	for e in events:
		type = e.type
		if type == "NoteOnEvent":
			pitch = e.pitch
			tick = e.type
			velocity = e.type
			channel = e.channel
			start = e.tick
			Notes[channel].append(pitch)
			Onsets[channel].append(start)
			#print start
		if type == "NoteOffEvent":
			#print "noteoff"
			pitch = e.pitch
			tick = e.type
			velocity = e.type
			channel = e.channel
			end = e.tick
			Durations[channel].append(end-start)
			if end-start<minduration:minduration = end-start
	return Notes,Durations,Onsets,minduration,resolution

def readMIDItype2(file):
	pattern = midi.read_midifile(file)
	events = pattern.trackpool

	#We assume a 4 channel MIDI file
	Notes = {0:[],1:[],2:[],3:[]}
	Durations = {0:[],1:[],2:[],3:[]}
	Onsets = {0:[],1:[],2:[],3:[]}
	minduration = 10000000
	ons = []
	for e in events:
		type = e.type
		if type == "NoteOnEvent":
			if e.velocity == 0: type = "NoteOffEvent"
		if type == "NoteOnEvent":
			pitch = e.pitch
			tick = e.type
			velocity = e.type
			channel = e.channel
			start = e.tick
			Notes[channel].append(pitch)
			Onsets[channel].append(start)
			ons.append(start)
			#print start
		if type == "NoteOffEvent":
			#print "noteoff"
			pitch = e.pitch
			tick = e.type
			velocity = e.type
			channel = e.channel
			end = e.tick
			start = ons[0]
			Durations[channel].append(end-start)
			ons.pop(0)
			if end-start<minduration:minduration = end-start
	return Notes,Durations,Onsets,minduration,_

def readTxtFile(file):
	f = open(file, 'r')
	text = f.read()
	f.close()
	return text

def filesInPath(mypath,ext = []):
	if len(ext) == 0:
		onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f!= ".DS_Store" ]
	else :
		onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) and f!= ".DS_Store" and f[-len(ext):] == ext]
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
	newS = []
	ids = []
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
		if type(seq) == ListType:
			f.write("".join(seq))
		if type(seq) == StringType or type(seq) == UnicodeType:
			f.write(seq)
		f.write("\n")
	f.close

def exportGroupOfSequencesToFASTAwithIDs(num_seq_MSA,ids,output):
	f = open(output, 'w')
	f.close
	f = open(output, 'a')
	for i,seq in enumerate(num_seq_MSA):
		f.write(">"+ids[i]+"\n")
		if type(seq) == ListType:
			f.write("".join(seq))
		if type(seq) == StringType or type(seq) == UnicodeType:
			f.write(seq)
		f.write("\n")
	f.close

