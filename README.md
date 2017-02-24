# BlastTools
## Purpose
A set of tools for preparing, running and processing BLAST searches 
## ExtractCDS
Extracts CDS sections of sequences.

Written in python3.  

Basic usage is:
```bash
ExtractCDS -i input.fas -o nucleotide
```
ExtractCDS takes two arguments (both of which are mandatory):

	-i (input fasta file)

	-o (output fasta format: nucleotide or protein)

NB: fasta file names MUST be in the form '>LabelX LabelY CDS=1-90', where the last space-delimited name element contains the coordinates of the CDS preceded by 'CDS=' (as output by gffread)
## OrthoFinder
Finds orthologues in a query set of proteins, based on similarity to a reference set of proteins and domain content.

Written in python3.

Basic usage is:
```bash
OrthoFinder -r reference.fasta -q query.fasta -g gene1 -d domain1 -o outputname
```
OrthoFinder takes five mandatory arguments:

-r (reference protein set in fasta format)

-q (query protein set in fasta format)

-g (gene of interest)

-d (domain to screen potential orthologues on)

-o (basename of output file)

OrthoFinder takes two optional arguments:

-c (number of threads to run analysis on, default=1)

--noblast (use existing blast databases to query a new gene of interest)

## ConservedOrthoFinder
Finds orthologues that are conserved as single copies across a number of protein sets, based on similarity to a reference set of proteins. Outputs separate fasta files for each single-copy orthologue found, and a fasta file (Concatenated.fasta) of all single-copy orthologues concatenated together (in the same order) for each protein set.

Written in python3.

Basic usage is:
```bash
ConservedOrthoFinder -r reference.fasta -i input1.fasta,input2.fasta -d directory1,directory2 -o output1,output2
```
ConservedOrthoFinder takes four mandatory arguments:

-r (reference protein set in fasta format)

-i (comma-separated list of query protein sets in fasta format)

-d (comma-separated list of working directories)

-o (comma-separated list of output filename stems)

ConservedOrthoFinder takes two optional arguments:

-c (number of threads to run analysis on, default=1)

-e (E-value for BLAST, default=0.001)

