#!/usr/bin/env python3

import argparse
import re
import sys
import os
import subprocess
import shutil
from Bio import SeqIO

#user argument parsing
parser=argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-d','--directory',type=str,help='Working directory')
parser.add_argument('-i','--input',type=str,help='Fasta file of all proteins for species of interest')
parser.add_argument('-r','--reference',type=str,help='Fasta file of all reference proteins to be found')
parser.add_argument('-e','--evalue',type=str,help='e-value threshold (default=0.1)')
parser.add_argument('-c','--cores',type=int,help='Number of cores to use (default=1)')
parser.add_argument('-o','--output',type=str,help='Output file name stem')
args=parser.parse_args()
#default arguments
Evalue=0.001
Cores=1
OutFile = "OutFile"
#working directory parsing
if args.directory is not None:
	DirectoryList = []
	for dir in args.directory.split(","):
		DirectoryList.append(dir)
		print("Working directory parsed: " + dir)
else:
	print("ERROR: working directory (-d) not specified")
	sys.exit(0)
#input fasta parsing
if args.input is not None:
	InputList = []
	for infile in args.input.split(","):
		InputList.append(infile)
		print("Input file parsed: " + infile)
else:
	print("ERROR: input fasta (-i) not specified")
	sys.exit(0)
#reference fasta parsing
if args.reference is not None:
	ReferenceFasta=args.reference
	print("Reference fasta = " + ReferenceFasta)
else:
	print("ERROR: reference fasta (-r) not specified")
	sys.exit(0)
#e-value parsing
if args.evalue is None:
	print("Using default e-value ("+str(Evalue)+")")
else:
	if float(args.evalue) > 0:
		Evalue=float(args.evalue)
		print("e-value set to "+str(Evalue))
	else:
		print("ERROR: e-value (-e) must be > 0")
#cores parsing
if args.cores is None:
	print("Using default number of cores ("+str(Cores)+")")
else:
	if int(args.cores) > 0:
		Cores=int(args.cores)
		print("Number of cores set to "+str(Cores))
	else:
		print("ERROR: cores (-c) must be > 0")
#output file name parsing
if args.output is not None:
	OutputList = []
	for outfile in args.output.split(","):
		OutputList.append(outfile)
		print("Output file name stored: "+OutFile)
else:
	print("ERROR: output file name (-o) not specified")
	sys.exit(0)

#log current directory (to change back to after all ConservedOrthoFinder runs finished, before running ConservedOrthoParser)
OriginalDir=os.getcwd()

def ConservedOrthoFinder(workdir,input,cores,reference,evalue,outfile):
	#change directory to working directory
	os.chdir(workdir)
	#make new directory to hold blast databases, outputs and results
	if os.path.isdir("ConservedOrthoFinder") is True:
		shutil.rmtree("ConservedOrthoFinder")
		print("Old ConservedOrthoFinder directory removed") 
	os.mkdir("ConservedOrthoFinder")
	#build DIAMOND blast database for input protein set
	cmd="diamond makedb --in " + input + " -d ./ConservedOrthoFinder/Input"
	subprocess.call(cmd,shell=True)
	print("BLAST database written for input")
	#blastp input using reference proteins
	cmd="diamond blastp -p " + str(cores) + " -d ./ConservedOrthoFinder/Input -q " + reference + " -o ./ConservedOrthoFinder/"+outfile+".out -f 6 qseqid evalue length sseqid -e " + str(evalue) + " --max-target-seqs 1 --more-sensitive"
	subprocess.call(cmd,shell=True)
	#collect top hit to each reference protein
	RefSeqs = []
	Hits = []
	for line in open('./ConservedOrthoFinder/'+outfile+'.out','r'):
		temp=line.split('\t')
		ref=temp[0]
		hit=temp[3].strip('\n')
		RefSeqs.append(ref)
		Hits.append(hit)
	print(str(len(Hits)) + ' hits parsed')
	#retrieve protein coding sequence for each hit
	Seqs = []
	for i in Hits:
		cmd='grep -A1 \"' + i + '\" ' + InputFasta + ' | grep -v \"' + i + '\"'
		seq = subprocess.getoutput(cmd)
		if seq is None:
			print('ERROR: sequence not found for '+i)
			sys.exit(0)
		elif seq=="/bin/sh:1:Syntaxerror:Unterminatedquotedstring":
			print('ERROR: sequence not found for '+i)
			sys.exit(0)
		else:	
			Seqs.append(seq.strip('\n').strip('*'))
	#check that the same number of refseqs, hits and sequences have been retrieved
	if not len(RefSeqs)==len(Hits)==len(Seqs):
		print('ERROR: different number of reference proteins (' + str(len(RefSeqs)) + ', hits (' + str(len(Hits)) + ') & sequences (' + str(len(Seqs)) + ')')
		sys.exit(0)
	#screen out any duplicate entries, making non-redundant (NR) lists of refseqs, hits and sequences
	RefSeqsNR = []
	HitsNR = []
	SeqsNR = []
	for i in range(len(RefSeqs)):
		if Seqs[i] not in SeqsNR:
			RefSeqsNR.append(RefSeqs[i])
			HitsNR.append(Hits[i])
			SeqsNR.append(Seqs[i])
	#check that the same number of NR refseqs, hits and sequences have been retrieved
	if not len(RefSeqsNR)==len(HitsNR)==len(SeqsNR):
		print('ERROR: different number of reference proteins (' + str(len(RefSeqsNR)) + ', hits (' + str(len(HitsNR)) + ') & sequences (' + str(len(SeqsNR)) + ')')
		sys.exit(0)
	#format output fasta file
	OutputFasta = ''
	for i in range(len(HitsNR)):
		templine = '>' + RefSeqsNR[i] + '_' + HitsNR[i] + '\n' + SeqsNR[i] + '\n'
		OutputFasta += templine
	#output fasta file of labelled orthologues
	outpath = os.getcwd()+'/ConservedOrthoFinder/'+OutFile+'.fasta'
	output = open(outpath,"wt")
	output.write(OutputFasta)
	output.close()
	print(str(len(HitsNR))+" orthologues written to "+outpath)
	return(outpath)

def ConservedOrthoParser(reference,inputs,inputnames):
	#change back to original directory
	os.chdir(OriginalDir)
	#split inputs and inputnames into lists for later iteration
	inputlist = []
	for i in inputs.split(","):
		if i != "":
			inputlist.append(i)
	inputnameslist = []
	for i in inputnames.split(","):
		if i != "":
			inputnameslist.append(i)
	#parse names of reference proteins
	RefNames = []
	for line in open(reference,"r"):
		if line.startswith(">"):
			RefNames.append(line.strip(">").strip("\n"))
	#go through each input file, generating a list of names of the sequences found in that input file, and adding them all to one central dict
	HitNames = {}
	HitSeqs = {}
	for i in range(len(inputnameslist)):
		Names = []
		Seqs = []
		for line in open(inputlist[i],"r"):
			if line.startswith(">"):
				Names.append(line.strip(">").strip("\n"))
			else:
				Seqs.append(line.strip("\n"))
		HitNames[inputnameslist[i]] = Names
		HitSeqs[inputnameslist[i]] = Seqs
	#go through the lists of sequence names for each input file, generating a count for each reference gene
	Presence = {}
	for inputname in HitNames:
		status = []
		for name in RefNames:
			found=0
			for hit in HitNames[inputname]:
				if re.search(name,hit):
					found+=1
			status.append(found)
			found=0
		Presence[inputname]=status
	#count up how many reference proteins are in each species, and how many species posses each reference protein
	RefTotals=[]
	for i in range(len(RefNames)):
		matchtotal=0
		for species in Presence:
			matchtotal+=Presence[species][i]
		RefTotals.append(matchtotal)
	SpeciesTotals=[]
	for species in inputnameslist:
		total=0
		for status in Presence[species]:
			total+=status
		SpeciesTotals.append(total)
	print("Total conserved orthologues found per species:")
	for i in range(len(inputnameslist)):
		print(inputnameslist[i]+":"+str(SpeciesTotals[i]))
	#assign reference proteins to either present in all, absent in at least one, or duplicated in at least one
	SingleCopy=[]
	AbsentInSome=[]
	Duplicated=[]
	for i in range(len(RefNames)):
		if RefTotals[i]==len(inputnameslist):
			SingleCopy.append(RefNames[i])
		elif RefTotals[i]<len(inputnameslist):
			AbsentInSome.append(RefNames[i])
		elif RefTotals[i]>len(inputnameslist):
			Duplicated.append(RefNames[i])
	print(str(len(SingleCopy))+" reference proteins present as single copies")
	print(str(len(AbsentInSome))+" reference proteins absent in at least one species")
	print(str(len(Duplicated))+" reference proteins duplicated in at least one species")
	#remove absent or duplicated sequences from the hit names and seqs dicts
	HitNamesSingle = {}
	HitSeqsSingle = {}
	for species in HitNames:
		singlenames=[]
		singleseqs=[]
		for i in range(len(HitNames[species])):
			for name in SingleCopy:
				if re.search(name,HitNames[species][i]):
					singlenames.append(HitNames[species][i])
					singleseqs.append(HitSeqs[species][i])
		HitNamesSingle[species]=singlenames
		HitSeqsSingle[species]=singleseqs
	#output fasta file for each single-copy sequence:
	for i in range(len(SingleCopy)):
		outfasta=''
		for j in HitNamesSingle:
			outfasta+=">"+HitNamesSingle[j][i]+"\n"+HitSeqsSingle[j][i]+"\n"
		output=open(SingleCopy[i]+".fasta","wt")
		output.write(outfasta)
		output.close()
	print("Individual fasta files written")
	#output concatenated fasta file of single-copy sequences for each species:
	HitConcatenation = {}
	for species in HitSeqsSingle:
		concatseq = ''
		for seq in HitSeqsSingle[species]:
			concatseq+=seq
		HitConcatenation[species]=concatseq
	outfasta=''
	for species in HitConcatenation:
		outfasta+=">"+species+"\n"+HitConcatenation[species]+"\n"
	output=open("Concatenated.fasta","wt")
	output.write(outfasta)
	output.close()
	print("Concatenated fasta written to Concatenated.fasta")
	return()

#find top hit to each reference protein in each input fasta
ResultsDirs=[]
for i in range(len(DirectoryList)):
	if not len(DirectoryList)==len(InputList)==len(OutputList):
		sys.exit(0)
	WorkDir = DirectoryList[i]
	print("Working directory changed: "+WorkDir)
	InputFasta = InputList[i]
	print("Input fasta changed: "+InputFasta)
	OutFile = OutputList[i] 
	print("Output file stem changed: "+OutFile)
	ResultsDirs.append(ConservedOrthoFinder(workdir=WorkDir,input=InputFasta,cores=Cores,reference=ReferenceFasta,evalue=Evalue,outfile=OutFile))

#parse results of ConservedOrthoFinder for each species
ResultsDirString = ''
for i in ResultsDirs:
	ResultsDirString+=i+","
InputNameString=''
for i in OutputList:
	InputNameString+=i+","
ConservedOrthoParser(reference=ReferenceFasta,inputs=ResultsDirString,inputnames=InputNameString)

