#!/usr/bin/env python
# coding: utf-8

# In[6]:


import ranking
import pandas as pd

datar = pd.read_csv("./Data/ProgramResults/OriginalPredictions/Ecoli/Final_sRNARFTarget_Ecoli.txt", header=0, delimiter = "\t") 
datar.iloc[:,-2] # Scores

rankr = list(ranking.Ranking(datar.iloc[:,-2], strategy=ranking.ORDINAL, start=1))
dfr = pd.DataFrame(rankr, columns = ['Rank', 'Score'])
dfr.to_csv('OrdinalRanked_Final_sRNARFTarget_Ecoli.txt', index = False, sep = '\t')


# In[ ]:




datac = pd.read_csv("./Data/ProgramResults/OriginalPredictions/Ecoli/Final_CopraRNA_Ecoli.txt", header=0, delimiter = "\t") 
datac.iloc[:,-2] # Scores

rankc = list(ranking.Ranking(datac.iloc[:,-2], strategy=ranking.ORDINAL, start=1))
dfc = pd.DataFrame(rankc, columns = ['Rank', 'Score'])
dfc.to_csv('OrdinalRanked_Final_CopraRNA_Ecoli.txt',index = False, sep = '\t')


# In[ ]:




datai = pd.read_csv("./Data/ProgramResults/OriginalPredictions/Ecoli/Final_IntaRNA_Ecoli.txt", header=0, delimiter = "\t") 
datai.iloc[:,-2] # Scores

ranki = list(ranking.Ranking(datai.iloc[:,-2], strategy=ranking.ORDINAL, start=1))
dfi = pd.DataFrame(ranki, columns = ['Rank', 'Score'])
dfi.to_csv('OrdinalRanked_Final_IntaRNA_Ecoli.txt',index = False, sep = '\t')

