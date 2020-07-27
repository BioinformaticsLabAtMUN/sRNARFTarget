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
  
    We created two python scripts for the understanding the predictions generated by sRNARFTarget. 
    sRNARFTarget_SHAP.py and sRNARFTarget_CP.py
    These are run after sRNARFTarget predictions are generated.
   
   #### 1. Instructions to run sRNARFTarget_SHAP
   
      1. Choose an sRNA-mRNA pair of interest from Prediction probabilities.csv file under sRNARFTargetResult folder.
      2. Run sRNARFTarget_SHAP using below command.
      
         python sRNARFTarget_SHAP.py sRNA_ID mRNA_ID
         
         Example usage: python sRNARFTarget_SHAP.py 'omrA' 'ompT'        
  
   #### 2. Instructions to run sRNARFTarget_CP
   
      1. For the same sRNA-mRNA pair that was chosen to run sRNARFTarget_SHAP, choose a feature/variable by looking at the plots generated by sRNARFTarget_SHAP.
      2. Run sRNARFTarget_CP using below command.
      
         python sRNARFTarget_CP.py sRNA_ID mRNA_ID feature_name
         
         Example usage: python sRNARFTarget_CP.py 'omrA' 'ompT' 'GCG'
         
   #### 3. Parameters
   
         1. 1st, 2nd and 3rd parameters are common for both programs.
         2. 1st parameter is python script name, sRNARFTarget_SHAP.py
         3. 2nd parameter is sRNA id same as in Prediction_probabilities.csv file under sRNARFTargetResult directory and should be in single quotes. For example: 'omrA'
         4. 3rd parameter is mRNA id same as in Prediction_probabilities.csv file under sRNARFTargetResult directory and should be in single quotes. For example: 'ompT'
         5. 4th parameter in sRNARFTarget_CP command is feature/variable name and should be in single quotes. For example : 'GCG'
      
   #### 4. Notes
         1. 2nd and 3rd parameter are ids of sRNA-mRNA pair taken from Prediction_probabilities.csv file under sRNARFTargetResult directory.
         2. sRNARFTarget_SHAP and sRNARFTarget_CP programs can be run for one sRNA-mRNA pair at a time.
