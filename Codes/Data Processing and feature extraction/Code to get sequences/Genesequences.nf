#!/usr/bin/env nextflow
#!/usr/bin/env edirect
//-------------------------channels----------------------------//
myFileChannel = Channel.fromPath('./foundsRNAsmRNAs.txt')
myFileChannel.into{channelv1; channelv2; channelv3}

channel2 = Channel.fromPath( './genseq.sh')

channel3 = Channel.fromPath( './DirgenomeSequences/*.fasta' )

channel4 = Channel.fromPath( './DirgenomeSequences/*.fasta' )

//-------------------------Process_1---------------------------//
process getAccessions{

  input:
    file p1f1 from channelv1

  output:
    file "*.txt" into p1result 

  script:
  """
  cut -f1 $p1f1 > accessions.txt
  """
}

p1result.into{P1A; P1B; P1W1}

//-------------------------Process_2---------------------------//
process getGenomeSequence{
 
 publishDir 'DirgenomeSequences'

  input:
    file p2f1 from P1A
    file p2f2 from channel2

  output:
    file "*.fasta" into p2result mode flatten

  script:
  """
  
  #cat $p2f1 | uniq > p2f3.txt
  awk '{print \$1}' $p2f1 | uniq > accession1.txt
  awk '{ if(\$1 != "Accession") print \$0;}' accession1.txt | sort | uniq > p2f3.txt
  n=1
  cat p2f3.txt | while read line; 
  do
    source ./$p2f2 
  n=\$((n+1))
  done 
  """  
}

p2result.collectFile().into{P2A; P2W1}

//-------------------------Process_3---------------------------//
process getGenomeLength{
 
  input:
    file p3f1 from P1B
    file p3f2 from P2A

  output:
    file "genomelength.txt" into p3result

  script:
  """
 IFS=\$'\n'; for next in \$(cat $p3f1); do esearch -db assembly -query "\$next [ACCN]" | elink -target nuccore -name assembly_nuccore_refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Slen Title; done >genomelength1.txt

 cat genomelength1.txt | cut -f1-2 > genomelength.txt
  """
}

p3result.into{P3A; P3B; P3W1}

//-------------------------Process_4---------------------------//
process getsRNADetails{
publishDir 'DirsRNADetails'
  input:
    file p4f1 from channelv2
    file p4f2 from P1W1
    file p4f3 from P2W1
    file p4f4 from P3W1

  output:
    file 'sRNAdetails.txt' into p4result

script:
"""
cut -f1,2 ${p4f1} > sRNAlist.txt

cat 'sRNAlist.txt' | while IFS=\$'\t' read a1 a2
do
esearch -db gene -query "\$a1[ACCN] AND \$a2[GENE]" < /dev/null | efetch -format docsum |
xtract -pattern DocumentSummary -element Name -block GenomicInfoType -element ChrAccVer ChrStart ChrStop > process1op1.txt
    cat 'process1op1.txt' | while read line
    do
    A="\$(cut -d\$'\t' -f1 <<<"\$line")"
        if [ \$A == \$a2 ];
        then
            echo \$line >> process1op2.txt
        else
            echo "\${line/\$A/\$a2}"  >> process1op2.txt
        fi
    done
done

awk '{if (\$3 > \$4) print \$2"\\t"\$4"\\t"\$3"\\t"\$1"\\t.\\t-"; else print \$2"\\t"\$3"\\t"\$4"\\t"\$1"\\t.\\t+"}' process1op2.txt > process1op3.txt

awk -F'\t' '{ print \$1"\\t"\$2"\\t"\$3"\\t"tolower(\$4)"\\t"\$5"\\t"\$6}' < process1op3.txt > process1op4.txt

awk 'FNR==NR{a[\$1,\$2]; next} (\$1,\$4) in a{print}' sRNAlist.txt process1op4.txt > sRNAdetails.txt

"""
}

p4result.set{P4A}

//-------------------------Process_5---------------------------//
process getsRNABedFile{
publishDir 'DirsRNABed'
container 'genomicpariscentre/bedtools'

  input:
  file p5f1 from P3A //genomelength
  file p5f2 from P4A //srnadetails

  output:
  file 'sRNAbedfile.bed' into p5result

  script:
  """
  bedtools slop -i ${p5f2} -g ${p5f1} -l 0 -r 1 > sRNAbedfile.bed
  """
}

p5result.set{P5A}

//-------------------------Process_6---------------------------//

process getsRNASequencesfrombed{
publishDir 'DirSequences'
container 'genomicpariscentre/bedtools'

  input:
    file p6f1 from P5A //bedfile
    file p6f2 from channel3 //genomefasta

  output:
    file 'sRNAsequences.txt' into p6result

  script:
  """
  bedtools getfasta -fi ${p6f2} -bed ${p6f1} -name -s -tab > genesequences1.fasta

  awk 'BEGIN {print "sRNA\tSequence"} {print}' genesequences1.fasta > sRNAsequences.txt
  """
}

p6result.collectFile(name: file("sRNAsequences.txt")).set{P6A}

//-------------------------Process_7---------------------------//
process getmRNADetails{
publishDir 'DirmRNADetails'
input:
file p7f1 from channelv3
file p7f2 from P6A

output:
file 'mRNAdetails.txt' into p7result

script:
"""
cut -f1,3 ${p7f1} > mRNAlist.txt

cat 'mRNAlist.txt' | while IFS=\$'\t' read a1 a2
do
esearch -db gene -query "\$a1[ACCN] AND \$a2[GENE]" < /dev/null | efetch -format docsum |
xtract -pattern DocumentSummary -element Name -block GenomicInfoType -element ChrAccVer ChrStart ChrStop > process2op1.txt
    cat 'process2op1.txt' | while read line
    do
    A="\$(cut -d\$'\t' -f1 <<<"\$line")"
        if [ \$A == \$a2 ];
        then
            echo \$line >> process2op2.txt
        else
            echo "\${line/\$A/\$a2}"  >> process2op2.txt
        fi
    done
done

awk '{if (\$3 > \$4) print \$2"\\t"\$4"\\t"\$3"\\t"\$1"\\t.\\t-"; else print \$2"\\t"\$3"\\t"\$4"\\t"\$1"\\t.\\t+"}' process2op2.txt > process2op3.txt

awk -F'\t' '{ print \$1"\\t"\$2"\\t"\$3"\\t"tolower(\$4)"\\t"\$5"\\t"\$6}' < process2op3.txt > process2op4.txt

awk 'FNR==NR{a[\$1,\$2]; next} (\$1,\$4) in a{print}' mRNAlist.txt process2op4.txt > mRNAdetails.txt

"""
}
p7result.set{P7A}
//-------------------------Process_8---------------------------//
process getmRNABedFile{
publishDir 'DirmRNABed'
container 'genomicpariscentre/bedtools'

input:
file p8f1 from P3B
file p8f2 from P7A

output:
file 'mRNAbedfile.bed' into p8result

script:
"""
bedtools slop -i ${p8f2} -g ${p8f1} -l 0 -r 1 > mRNAbedfile.bed
"""
}

p8result.set{P8A}

//-------------------------Process_9---------------------------//

process getmRNASequencesfrombed{

container 'genomicpariscentre/bedtools'

input:
file p9f1 from P8A
file p9f2 from channel4

output:
file 'mRNAsequences.txt' into p9result

script:
"""
bedtools getfasta -fi ${p9f2} -bed ${p9f1} -name -s -tab > genesequences2.fasta

awk 'BEGIN {print "mRNA\tSequence"} {print}' genesequences2.fasta > mRNAsequences.txt
"""
}

p9result.collectFile(name: file("mRNAsequences.txt"))
