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
parser.add_argument('-r','--reference',type=str,help='Fasta file of all proteins for reference genome')
parser.add_argument('-q','--query',type=str,help='Fasta file of all proteins for query genome')
parser.add_argument('-e','--evalue',type=str,help='e-value threshold (default=0.1)')
parser.add_argument('-g','--gene',type=str,help='Gene of interest (case insensitive)')
parser.add_argument('-c','--cores',type=int,help='Number of cores to use (default=1)')
parser.add_argument('-o','--outfile',type=str,help='Output file name')
parser.add_argument('-d','--domains',type=str,help='Domains to screen hits on (only hits containing all domains will be retained)')
parser.add_argument('--noblast',help="Use existing all vs all blast databases",action="store_true")
args=parser.parse_args()
#default arguments
Evalue=0.1
Cores=1
Domains=[]
NoBLAST=False
#reference fasta parsing
if args.reference is not None:
	ReferenceFasta=args.reference
	print("Reference fasta = " + ReferenceFasta)
else:
	print("ERROR: reference fasta (-r) not specified")
	sys.exit(0)
#query fasta parsing
if args.query is not None:
	QueryFasta=args.query
	print("Query fasta = " + QueryFasta)
else:
	print("ERROR: query fasta (-q) not specified")
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
#genes of interest parsing
if args.gene is not None:
	print("Gene of interest = " + args.gene)
	Gene=args.gene
else:
	print("ERROR: gene of interest (-g) not specified")
	sys.exit(0)
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
if args.outfile is not None:
	OutFile=args.outfile
	print("Output file name = " + OutFile)
else:
	print("ERROR: output file name (-o) not specified")
	sys.exit(0)
#domain parsing
if args.domains is not None:
	print("Domains to screen on:")
	for i in range(len(args.domains.split(','))):
		Domains.append(args.domains.split(',')[i])
		print(args.domains.split(',')[i])
else:
	print("No domain screening specified")
#noblast parsing
if args.noblast:
	NoBLAST=True
	print("Using existing BLAST databases")
else:
	print("Creating new BLAST databases")

if NoBLAST==False:
	#make new directory to hold blast databases, outputs and results
	if os.path.isdir("OrthoFinder") is True:
		shutil.rmtree("OrthoFinder")
		print("Old OrthoFinder directory removed") 
	os.mkdir("OrthoFinder")
	#build DIAMOND blast databases for query & reference protein sets
	cmd="diamond makedb --in " + QueryFasta + " -d ./OrthoFinder/Query"
	subprocess.call(cmd,shell=True)
	cmd="diamond makedb --in " + ReferenceFasta + " -d ./OrthoFinder/Reference"
	subprocess.call(cmd,shell=True)
	print("BLAST databases written for query and reference")
	#blastp reference using query
	cmd="diamond blastp -p " + str(Cores) + " -d ./OrthoFinder/Reference -q " + QueryFasta + " -o ./OrthoFinder/QueryOnReference.out -f 6 qseqid evalue length sseqid -e " + str(Evalue) + " --max-target-seqs 1 --more-sensitive"
	subprocess.call(cmd,shell=True)
	#blastp query using reference
	cmd="diamond blastp -p " + str(Cores) + " -d ./OrthoFinder/Query -q " + ReferenceFasta + " -o ./OrthoFinder/ReferenceOnQuery.out -f 6 qseqid evalue length sseqid -e " + str(Evalue) + " --max-target-seqs 1 --more-sensitive"
	subprocess.call(cmd,shell=True)
#collect top hits to reference for each query sequence
QueryOnRef = {}
for line in open('./OrthoFinder/QueryOnReference.out','r'):
	temp=line.split('\t')
	QueryName=temp[0]
	RefName=temp[3].strip('\n')
	QueryOnRef[QueryName]=RefName
print(str(len(QueryOnRef)) + ' Query-Reference pairs parsed')
#find match to gene of interest in reference hits
QueryOnRefHits = []
QueryOnRefOrtho = []
for i in QueryOnRef:
	if re.search(Gene,QueryOnRef[i],flags=re.IGNORECASE):
		QueryOnRefHits.append(i)
		QueryOnRefOrtho.append(Gene)
#collect top hits to query for each reference sequence
RefOnQuery = {}
for line in open('./OrthoFinder/ReferenceOnQuery.out','r'):
	temp=line.split('\t')
	RefName=temp[0]
	QueryName=temp[3].strip('\n')
	RefOnQuery[QueryName]=RefName
print(str(len(RefOnQuery)) + ' Reference-Query pairs parsed')
#find match to gene of interest in reference hits
RefOnQueryHits = []
RefOnQueryOrtho = []
for i in RefOnQuery:
	if re.search(Gene,RefOnQuery[i],flags=re.IGNORECASE):
		RefOnQueryHits.append(i)
		RefOnQueryOrtho.append(Gene)	
#combine QueryOnRefHits & RefOnQueryHits into one non-redundant list 
nrHits = []
nrOrtho = []
for i in range(len(QueryOnRefHits)):
	if QueryOnRefHits[i] not in nrHits:
		nrHits.append(QueryOnRefHits[i])
		nrOrtho.append(QueryOnRefOrtho[i])
for i in range(len(RefOnQueryHits)):
	if RefOnQueryHits[i] not in nrHits:
		nrHits.append(RefOnQueryHits[i])
		nrOrtho.append(RefOnQueryOrtho[i])
#retrieve protein coding sequences for the nrHits
nrSeqs = []
for hit in nrHits:
	found=0
	cmd='grep -A1 \'' + hit + '\' ' + QueryFasta + ' | grep -v \'' + hit + '\''
	seq = subprocess.getoutput(cmd)
	if seq is not None:
		found=1
		nrSeqs.append(seq.strip('\n').strip('*'))
	if found==0:
		print('ERROR: sequence not found for hit ' + hit)

#check that the same number of names and sequences have been retrieved before output
if not len(nrHits)==len(nrOrtho)==len(nrSeqs):
	print('ERROR: different number of hits (' + str(len(nrHits)) + '), orthologues (' + str(len(nrOrtho)) + ') & sequences (' + str(len(nrSeqs)) + ')')
	sys.exit(0)
#format output fasta file
Ortho=[]
for i in range(len(nrHits)):
	Ortho.append(nrOrtho[i]+'_'+nrHits[i])
print('QueryOnRefHits (' + str(len(QueryOnRefHits)) + ' entries) & RefOnQueryHits (' + str(len(RefOnQueryHits)) + ' entries) combined into Ortho (' + str(len(Ortho)) + ' entries)')
OrthoFasta = ''
for i in range(len(nrHits)):
	templine = '>' + Ortho[i] + '\n' + str(nrSeqs[i]) + '\n'
	OrthoFasta += templine
#output fasta file of labelled orthologues
output = open('./OrthoFinder/'+OutFile+'.fasta',"wt")
output.write(OrthoFasta)
output.close()
if len(Domains)==0:
	print("Orthologues written to ./OrthoFinder/" + OutFile + ".fasta")
else:
	#screen out hits without domains
	#identify domains in each hit using InterProScan with Pfam database
	cmd = "interproscan.sh -appl Pfam -i ./OrthoFinder/" + OutFile + ".fasta -o ./OrthoFinder/" + OutFile + ".fasta.tsv -f tsv"
	subprocess.call(cmd,shell=True)
	#collect list of domains for each hit
	OrthoDomains = []
	for hit in Ortho:
		tempdomains=[]
		for line in open("./OrthoFinder/" + OutFile + ".fasta.tsv"):
			if hit in line:
				tempdomains.append(line.split('\t')[5])
		OrthoDomains.append(tempdomains)	
	#identify which hits have the domains of interest
	OrthoVerified = []
	OrthoVerifiedSeqs = []
	for i in range(len(Ortho)):
		#count of domains found in this hit
		match=0
		for orthodomain in OrthoDomains[i]:
			for domain in Domains:
				if re.search(domain,orthodomain,flags=re.IGNORECASE):	
					match+=1
		if match>0:
			OrthoVerified.append(Ortho[i])
			OrthoVerifiedSeqs.append(nrSeqs[i])
			print(Ortho[i] + " added to verified list")
	#output fasta of domain-verified hits
	OrthoVerifiedFasta = ''
	for i in range(len(OrthoVerified)):
		templine = '>' + OrthoVerified[i] + '\n' + str(OrthoVerifiedSeqs[i]) + '\n'
	OrthoVerifiedFasta += templine
	output = open("./OrthoFinder/" + OutFile + '_verified.fasta',"wt")
	output.write(OrthoVerifiedFasta)
	output.close()
	print("Verified orthologues written to ./OrthoFinder/"+OutFile+"_verified.fasta")

