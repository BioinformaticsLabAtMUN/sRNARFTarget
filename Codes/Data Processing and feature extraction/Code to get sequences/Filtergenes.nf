#!/usr/bin/env nextflow

//-------------------------channel---------------------------//

FileChannel = Channel.fromPath( './Data.txt').into{channel1; channel2; channel3}
//-------------------------Process_1---------------------------//
process checkifsRNAsFound{

input:
file p1f1 from channel1

output:
file 'sRNAcount.txt' into p1result

script:
"""
cat ${p1f1} | cut -f1,2 | while IFS=\$'\t' read f1 f2
do
	esearch -db gene -query "\$f1[ACCN] AND \$f2[GENE]" < /dev/null| xtract -pattern ENTREZ_DIRECT -element Count 
done > sRNAcount.txt

"""
}

p1result.into{P1A; P1B}

//-------------------------Process_2---------------------------//
process checkifmRNAsFound{

input:
file p2f1 from channel2
file p2f2 from P1A

output:
file 'mRNAcount.txt' into p2result

script:
"""
cat ${p2f1} | cut -f1,3 | while IFS=\$'\t' read f1 f3
do

	esearch -db gene -query "\$f1[ACCN] AND \$f3[GENE]" < /dev/null| xtract -pattern ENTREZ_DIRECT -element Count 
done > mRNAcount.txt

"""
}

p2result.into{P2A; P2W1}

//-------------------------Process_3---------------------------//
process getfilteredsRNAmRNA{

publishDir 'DirfilteredgenesStarothers'

input:
file p3f1 from channel3
file p3f2 from P1B	   
file p3f3 from P2A	   

output:
file '*.txt' into p3result mode flatten

script:
"""
paste ${p3f1} ${p3f2} ${p3f3}| awk '{print \$1,\$2,\$3,\$4,\$5}' > process3op1.txt

awk '(\$4 != 0) && (\$5 != 0 )' process3op1.txt > sRNAmRNAwithcount1.txt

tr ' ' '\t' < sRNAmRNAwithcount1.txt | cut -f1-3 > foundsRNAsmRNAs.txt

awk '(\$4 == 0) || (\$5 == 0 )' process3op1.txt > sRNAmRNAwithcount0.txt

tr ' ' '\t' < sRNAmRNAwithcount0.txt | cut -f1-3 > notfoundsRNAsmRNA.txt
"""
}

p3result.collectFile(name: file("foundsRNAsmRNAs.txt")
