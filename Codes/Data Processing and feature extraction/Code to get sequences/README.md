# Filtergenes.nf

## Introduction

This pipeline filter the file containing sRNA-mRNA pairs and generates an output file containing the found sRNA-mRNA pairs in the NCBI pubmed. First it checks for sRNA counts and then for mRNA using the accession number and gene name. Finally, it creates the output file called *foundsRNAsmRNAs.txt* using the sRNA and mRNA count files.

## Input File:

The input file contains the accession numbers, name of sRNAs, and name of mRNAs in lower case.
Example:
```
Accession     sRNA  mRNA
NC_017501.1   sdsr  glgX
```
# Genesequences.nf

## Introduction

This pipeline gets the sequence of sRNA/mRNA found in the Filtergenes pipeline. It has in total 9 process to get the sequences. Process 1 outputs the file containing the accession number from the input file. Process 2 and 3 gets the sequence of genome and genome length, respectively. From process 4 to 6 it gets the sequences of all the sRNAs and from 7 to 9 it gets the sequences of all the mRNAs. Process 6 outputs the text file containing the sRNA names with their sequences and process 9 outputs the mRNA names with their sequences.

## Input File:

The output file *foundsRNAsmRNAs.txt* generated in Filtergenes.nf that contains the accession numbers, name of found sRNAs, and name of found mRNAs.