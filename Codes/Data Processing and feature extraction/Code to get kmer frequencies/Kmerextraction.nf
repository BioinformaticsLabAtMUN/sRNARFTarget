#!/usr/bin/env nextflow

//Usage: nextflow run tutorial.nf --k '3'

Channel1 = Channel.fromPath( './sRNASequences.txt' )
Channel2 = Channel.fromPath( './mRNASequences.txt' )

params.k = 0

process getsRNAFeatures{

  input:
    file File6 from Channel1

  output:
    file 'sRNA_kmer.txt' into process5result

script:
"""
#!/usr/bin/env python3
  
from pprint import pprint
from skbio import Sequence
from itertools import product
from sklearn.utils import shuffle
import pandas as pd
import numpy as np


data1 = pd.read_csv('$File6', sep='\t')
df = pd.DataFrame(data = data1)

name = df.iloc[:, :-1].values #All columns but last
sequences = df.iloc[:,  -1].values #Only last column

seqarr = []

for item in sequences:
  number = 1
  s = Sequence(item)
  freqs = s.kmer_frequencies(${params.k}, relative=True, overlap=True)
  seqarr.append(freqs)
  number = number + 1

def all_kmer_subsets(ss=["A", "T", "G", "C"]):
  return [''.join(p) for p in product(ss, repeat=${params.k})]

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

df3.to_csv('sRNA_kmer.txt', header=True, index=False, sep='\t', mode='a')

"""
}

//Collect file
process5result
.collectFile(name: file("sRNA${params.k}_merfeatures.csv"))
.into{setResult5;setResult55}

process getmRNAFeatures{


  input:
    file File10 from Channel2

  output:
    file 'mRNA_kmer.txt' into process9result

script:
"""
#!/usr/bin/env python3
  
from pprint import pprint
from skbio import Sequence
from itertools import product
import pandas as pd
import numpy as np

data1 = pd.read_csv('$File10', sep='\t', encoding = 'unicode_escape')
df = pd.DataFrame(data = data1)

name = df.iloc[:, :-1].values #All columns but last
sequences = df.iloc[:,  -1].values #Only last column

seqarr = []

for item in sequences:
  number = 1
  s = Sequence(item)
  freqs = s.kmer_frequencies(${params.k}, relative=True, overlap=True)
  seqarr.append(freqs)
  number = number + 1

def all_kmer_subsets(ss=["A", "T", "G", "C"]):
  return [''.join(p) for p in product(ss, repeat=${params.k})]

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

df3.to_csv('mRNA_kmer.txt', header=True, index=False, sep='\t', mode='a')

"""
  
}

//Collect file
process9result
.collectFile(name: file("mRNA${params.k}_merfeatures.csv"))
.into{setResult9;setResult99}



process mergesRNAmRNAfeaturefiles{

  input:
    file sRNAfile from setResult5
    file mRNAfile from setResult9

  output:
    file 'process10op3mer.csv' into process10result

  script:
  """
#!/usr/bin/env Rscript
  rm(list=ls());

  sRNAFile <- "$sRNAfile"
  mRNAFile <- "$mRNAfile"

  sRNAs <- read.table(sRNAFile, sep="\t", header=TRUE)
  mRNAs <- read.table(mRNAFile, sep="\t", header=TRUE)
  finalResult <- cbind(sRNAs, mRNAs)
  write.table(finalResult,file = "process10op3mer.csv", sep="\t",quote = FALSE, col.names=TRUE, row.names = FALSE);
  """
  
}

//Collect file
process10result
.collectFile(name: file("${params.k}_merFinalFeatures.csv"))
.set{setResult10}

process getKmerdifference{

  input:
    file skmerfile from setResult55
    file mkmerfile from setResult99

  output:
    file 'kmerdifference.txt' into process11result

  script:
  """
  #!/usr/bin/env python3

  import pandas as pd
  sRNA = pd.read_csv('$skmerfile', sep='\t')
  mRNA = pd.read_csv('$mkmerfile', sep='\t')
  sRNA
  mRNA
  #output8 = mRNA - sRNA
  output8 = mRNA.subtract(sRNA)
  output8.to_csv('kmerdifference.txt', header=False, index=False, sep='\t', mode='a')
  """
  
}

//Collect file
process11result
.collectFile(name: file("TrinucleotidesDifference.csv"))
.set{setResult11}
