#!/bin/bash

cat 'p2f3.txt' | while read line

do 
   esearch -db nucleotide -query "$line [ACCN]" < /dev/null | efetch -format fasta > "$line.fasta"
done 



