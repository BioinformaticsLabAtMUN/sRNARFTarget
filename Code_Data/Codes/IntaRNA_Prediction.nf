#!/usr/bin/env nextflow

//-------------------------channels----------------------------//


channelsrna = Channel.fromPath('./1srna.fasta')
channelmrna = Channel.fromPath('./Multocida1804mRNAsequences.fasta')
channelrscript = Channel.fromPath('./IntaRNA_CSV_p-value.R')


//------------------------Process_1---------------------------//

process getIntaRNAPredictions{

input:
file sfile from channelsrna
file mfile from channelmrna

output:
file '*.csv' into process1result

script:
"""
IntaRNA -q $sfile -t $mfile --outMode=C --outCsvCols=id2,id1,E,P_E --out=PCC_IntaRNA_Result_New.csv

"""
}
process1result.collectFile(name: file("PCC_IntaRNA_Result_New.csv")).set{setResult1}

//-------------------------Process_2---------------------------//

process getIntaRNAP_values{

input:
file efile from setResult1
file rfile from channelrscript

output:
file '*.csv' into process2result

script:
"""
Rscript --vanilla $rfile $efile PCC_IntaRNA_Result_pvalue_New.csv E

"""
}

//Collect file
process2result.collectFile(name: file("Multocida_IntaRNA_Result_pvalue_locus_tag.csv")).set{setResult2}



workflow.onComplete {
println(
"""
Pipeline execution summary
---------------------------
Run as : ${workflow.commandLine}
Completed at: ${workflow.complete}
Duration : ${workflow.duration}
Success : ${workflow.success}
workDir : ${workflow.workDir}
exit status : ${workflow.exitStatus}
""")
}

