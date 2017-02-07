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

OrthoFinder takes one optional argument:

-c (number of threads to run analysis on, default=1)
