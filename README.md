# sRNARFTarget: A machine learning-based approach for sRNA Target Prediction #
  
  ## Introduction

    This repository contains all the codes, data, results and supplementary files related to the 'sRNARFTarget' program for sRNA target prediction.
  
  ## Instructions to run sRNARFTarget
  
    1. Clone the repository.
    2. Keep the sRNA and mRNA fasta files of interest in the sRNARFTarget-master folder.
    3. Set 'sRNARFTarget-master' as the current working directory.
    4. Run the below command to run sRNARFTarget. Replace sRNA.fasta (--s parameter) and mRNA.fasta (--m parameter) with the desired fasta file of sRNAs and mRNA present in the sRNARFTarget-master directory.
   
  ## Command to run sRNARFTarget
    nextflow run sRNARFTarget.nf --s sRNA.fasta --m mRNA.fasta
   
  ## Creation of all possible sRNA-mRNA pairs
  
    sRNARFTarget creates all possible pairs from the input sRNA and mRNA sequences. Each sRNA is paired with all mRNAs. 
    For example, if the input sRNA file has 5 sRNA sequences and mRNA file has 9 mRNA sequences, then it will create 45 sRNA-mRNA pairs, 9 pairs for each sRNA.
 
  ## sRNARFTarget Results
    
    1. When the program has finished execution, it creates a directory 'sRNARFTargetResult' with two files.
    2. Prediction_probabilities.csv: this file is the sRNARFTarget result file and contains results sorted by predicted interaction probability from high to low, rounded to five decimals. It contains three columns, sRNA_ID, mRNA_ID and Prediction_Probability.
    3. FeatureFile.csv: this file contains features for all the sRNA-mRNA pairs. This file consists of 66 columns. The first two columns are sRNA_ID and mRNA_ID.
       The remaining 64 columns are corresponding trinucleotide frequency difference of sRNA-mRNA pairs. This file is later used by sRNARFTarget interpretability scripts.

  ## sRNARFTarget Predictions Interpretation
