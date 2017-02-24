#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO

#user argument parsing
parser=argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i','--infile',type=str,help='Input fasta file')
parser.add_argument('-o','--outputformat',type=str,help='Specify output fasta as either nucleotide (\'nucleotide\') or protein (\'protein\')')
args=parser.parse_args()
#fasta input file parsing
FastaFile=args.infile
if FastaFile is not None:
	print("Input fasta = " + FastaFile)
else:
	print("ERROR: no input fasta (-i) specified")
	sys.exit(0)
#output format parsing
OutFormat='nucleotide'
if args.outputformat != 'nucleotide':
	if args.outputformat=='protein':
		OutFormat='protein'
		print('Output format = protein')
	else:
		print('ERROR: invalid option specified for output format (-i) - valid options are \'nucleotide\' or \'protein\'')
		sys.exit(0)
else:
	print('Output format = nucleotide')
#read each name
Names=[]
StartPos=[]
EndPos=[]
for line in open(FastaFile):
	if line.startswith('>'):
		temp=line.split(' ')
		tempname = ''
		if temp[-1].startswith('CDS='):
			for i in range(len(temp)-1):
				tempname+=temp[i].lstrip('>')+'_'
			Names.append(tempname.strip('_'))
			positions=(temp[-1].split('-'))
			StartPos.append(int(positions[0].strip('CDS=')))
			EndPos.append(int(positions[1]))	
		else:
			print('ERROR: no CDS coordinates found in name ' + line)
			sys.exit(0)	
#read each sequence
Seqs=[]
for record in SeqIO.parse(FastaFile,'fasta'):
	Seqs.append(record.seq)
if len(Names)!=len(Seqs):
	print("ERROR: not all sequences have a name (" + str(len(Names)) + " names & " + str(len(Seqs)) + " sequences)")
	sys.exit(0)
#extract CDS and translation from each seq
CDS=[]
Proteins=[]
for i in range(len(Seqs)):
	CDS.append(Seqs[i][StartPos[i]-1:EndPos[i]])
	Proteins.append(Seqs[i][StartPos[i]-1:EndPos[i]].translate())
#find genes with stop codons in their CDS translation
StopSeqs=[]
for i in range(len(Names)):
	#if there is a stop codon anywhere in the protein...
	if Proteins[i].find('*')!=-1:
		#strip off the last aa (which could be a genuine stop codon) and retest protein for stops
		temp=Proteins[i].strip('*')
		#if stops still present (internal stops), record the name of that sequence
		if temp.find('*')!=-1:
			StopSeqs.append(Names[i])
if len(StopSeqs)>0:
	print('Stop codon found in CDS of ' + str(len(StopSeqs)) + '/' + str(len(Names)) + ' sequences')
#write cleaned-up CDS to file
CDSfasta = ''
if OutFormat=='nucleotide':
	for i in range(len(Names)):
		if Names[i] not in StopSeqs:
			CDSfasta += '>' + Names[i] + '\n' + str(CDS[i]) + '\n'
	if FastaFile.endswith('.fas'):
		OutName=FastaFile.replace('.fas','_clean_nucl.fas')
	elif FastaFile.endswith('.fasta'):
		OutName=FastaFile.replace('.fasta','_clean_nucl.fasta')
	else:
		OutName='CDS_clean_nucl.fasta'
if OutFormat=='protein':
	for i in range(len(Names)):
		if Names[i] not in StopSeqs:
			CDSfasta += '>' + Names[i] + '\n' + str(Proteins[i]) + '\n'
	if FastaFile.endswith('.fas'):
		OutName=FastaFile.replace('.fas','_clean_protein.fas')
	elif FastaFile.endswith('.fasta'):
		OutName=FastaFile.replace('.fasta','_clean_protein.fasta')
	else:
		OutName='CDS_clean_nucl.fasta'
OutFile=open(OutName,'wt')
OutFile.write(CDSfasta)
OutFile.close()
print('Cleaned CDS file written to ' + OutName)

