#!/usr/bin/env nextflow

//-------------------------channels----------------------------//
params.s = "srna.fasta"
params.m = "mrna.fasta"

srna = file("${params.s}")
mrna = file("${params.m}")

ldedmodel = Channel.fromPath('./PickledModelData/RFModel/sRNARFTargetModel.pickle')
if(!srna.exists()) {
exit 1, "The specified input sRNA fasta file does not exist: ${params.s}"
}

if(!mrna.exists()) {
exit 1, "The specified input mRNA fasta file does not exist: ${params.m}"
}

//------------------------Process_1---------------------------//

process createAllPossiblePairs{

input:
file afile from srna
file bfile from mrna

output:
file 'pairs_names_seqs.txt' into process1result mode flatten

script:
"""
#!/usr/bin/env python3
from Bio import SeqIO

ns = open('pairs_names_seqs.txt', 'w')

ns.write("sRNA" + "\t"+ "mRNA"+ "\t" + "sRNA_Sequence" + "\t" + "mRNA_Sequence")
ns.write('\\n')

for record1 in SeqIO.parse("$srna", "fasta"):
    for record2 in SeqIO.parse("$mrna", "fasta"):
        ns.write(record1.id + "\t"+ record2.id + "\t" + str(record1.seq) + "\t" + str(record2.seq))
        ns.write('\\n')
ns.close()
"""
}
process1result.collectFile(name: file("pairs_names_seqs.txt")).into{setResult1; setResult11; setResult111}
//-------------------------Process_2---------------------------//

process getsRNATrinucleotidesFrequncies{

input:
file nsfile from setResult1

output:
file 'sRNA_3mer.txt' into process2result

script:
"""
#!/usr/bin/env python3
from pprint import pprint
from skbio import Sequence
from itertools import product
from sklearn.utils import shuffle
import pandas as pd
import numpy as np


data1 = pd.read_csv('$nsfile', sep='\t', header=0)
df = pd.DataFrame(data = data1)
sequences = df.iloc[:, 2].values #third column sRNA sequences

seqarr = []

for item in sequences:
    number = 1
    s = Sequence(item)
    freqs = s.kmer_frequencies(3, relative=True, overlap=True)
    seqarr.append(freqs)
    number = number + 1

def all_kmer_subsets(ss=["A", "T", "G", "C"]):
    return [''.join(p) for p in product(ss, repeat=3)]

kmer_combinations = all_kmer_subsets()

df1 = pd.DataFrame(seqarr) #convert dicionary to dataframe
rows = len(df1.index)
cols = len(kmer_combinations)
d = np.zeros(shape=(rows,cols))
df2 = pd.DataFrame(data = d, columns=kmer_combinations)
df3 = pd.DataFrame()

for col in kmer_combinations:
    if col in df1.columns:
        df3[col] = df1[col]
    else:
        df3[col] = df2[col]

df3 = df3.fillna(0) #fill empty columns with zero (replace NaN with 0)

df3.to_csv('sRNA_3mer.txt', header=True, index=False, sep='\t', mode='a')

"""
}

//Collect file
process2result.set{setResult2}
//-------------------------Process_2---------------------------//
process getmRNATrinucleotidesFrequncies{

input:
file sRNas from setResult2
file ns2file from setResult11

output:
file '3merdifference.txt' into process3result

script:
"""
#!/usr/bin/env python3

from pprint import pprint
from skbio import Sequence
from itertools import product
import pandas as pd
import numpy as np

data1 = pd.read_csv('$ns2file', sep='\t', header=0)
df = pd.DataFrame(data = data1)

sequences = df.iloc[:, 3].values #fourth column mRNA sequences

seqarr = []

for item in sequences:
    number = 1
    s = Sequence(item)
    freqs = s.kmer_frequencies(3, relative=True, overlap=True)
    seqarr.append(freqs)
    number = number + 1

def all_kmer_subsets(ss=["A", "T", "G", "C"]):
    return [''.join(p) for p in product(ss, repeat=3)]

kmer_combinations = all_kmer_subsets()

df1 = pd.DataFrame(seqarr) #convert dicionary to dataframe
rows = len(df1.index)
cols = len(kmer_combinations)
d = np.zeros(shape=(rows,cols))
df2 = pd.DataFrame(data = d, columns=kmer_combinations)
df3 = pd.DataFrame()

for col in kmer_combinations:
    if col in df1.columns:
        df3[col] = df1[col]
    else:
        df3[col] = df2[col]

df3 = df3.fillna(0) #fill empty columns with zero (replace NaN with 0)

df3.to_csv('mRNA_3mer.txt', header=True, index=False, sep='\t', mode='a')

sRNA = pd.read_csv('$sRNas', sep='\t')
sRNAdf  = pd.DataFrame(data = sRNA)

mRNA = pd.read_csv('mRNA_3mer.txt', sep='\t')
mRNAdf  = pd.DataFrame(data = mRNA)

output8 = mRNAdf.subtract(sRNAdf)
output8.to_csv('3merdifference.txt', header=True, index=False, sep='\t', mode='a')
"""

}

//Collect file
process3result.into{setResult3; setResult33}

//-------------------------Process_3---------------------------//

process runRandomForestModel{

input:
file newdata from setResult3
file lrf from ldedmodel

output:
file 'Results_pred_probs.txt' into process5result

script:
"""
#!/usr/bin/env python3
import pandas as pd
import numpy  as np
import pickle

def pred_prob(testdf):
    # load the saved random forest model from disk
    loaded_RFmodel = pickle.load(open('$lrf', 'rb'))
    #predict probabilities for class 1
    predict_proba = loaded_RFmodel.predict_proba(testdf)
    print(predict_proba)
    #write probabilities to file
    for i in predict_proba:
        with open('Results_pred_probs.txt','a') as fd:
            fd.write(str(i[1])+"\\n")
    print(predict_proba[:, 1])
    return predict_proba[:, 1]

testdata = pd.read_csv('$newdata', sep='\t', header=0)
testdf = pd.DataFrame(data = testdata)
testdf = testdf.fillna(0)

pred_prob(testdf)
"""
}
process5result.set{setResult5}

//-------------------------Process_5---------------------------//

process generateSortedResultFile{
publishDir 'sRNARFTargetResult'
input:
file mlfile from setResult5
file ns3file from setResult111
file difffile from setResult33

output:
file '*.csv' into process6result

script:
"""
#!/usr/bin/env python3
import pandas as pd

#Generate sorted prediction result file
df1 = pd.read_csv('$ns3file', sep='\t', header=0)
df2 = pd.read_csv('$mlfile', sep='\t', header=None)
df3 = pd.DataFrame(data=df1.iloc[:, 0:2].values,columns=['sRNA', 'mRNA']).assign(Prediction_Probability=df2.round(5))

df4 = df3.sort_values('Prediction_Probability',ascending=False)
df4.to_csv('Prediction_probabilities.csv', sep='\t', index=False)

#Generate feature file with pair ids to be used for interpretability later
dfp61 = pd.read_csv('$difffile', sep='\t',header = 0) 
dfp63 = pd.concat([df1.iloc[:, 0:2], dfp61], axis = 1)
dfp63.to_csv('FeatureFile.csv', header = True, sep='\t', index=False)

"""
}



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
