#!/usr/bin/env python

#take 10% top predictions from all not srna wise and count the True positives, calculate % of tp and plot

import pandas as pd
ecolirf = pd.read_csv('./sRNARFTarget_Ecoli.txt', header = 0,sep = '\t')
ecolirf.head()

ecolirfrows = (ecolirf.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
erfdftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((ecolirfrows*i)/100)
    nows =  ecolirf.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/101
    erfdftoplot = erfdftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
erfdftoplot
erfdftoplot = erfdftoplot.round(2)
print(erfdftoplot)
Ecolirftolineplot = erfdftoplot[['%ofpredictions','%ofTruePositives']]
Ecolirftolineplot



import pandas as pd
ecolicopra = pd.read_csv('./CopraRNA_Ecoli.txt', header = 0,sep = '\t')
ecolicopra.head()

ecolicoprarows = (ecolicopra.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
ecopradftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((ecolicoprarows*i)/100)
    nows =  ecolicopra.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/101
    ecopradftoplot = ecopradftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
ecopradftoplot
ecopradftoplot = ecopradftoplot.round(2)
print(ecopradftoplot)
Ecolicopratolineplot = ecopradftoplot[['%ofpredictions','%ofTruePositives']]
Ecolicopratolineplot

import pandas as pd
ecoliinta = pd.read_csv('./IntaRNA_Ecoli.txt', header = 0,sep = '\t')
ecoliinta.head()

ecoliintarows = (ecoliinta.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
eintadftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((ecoliintarows*i)/100)
    nows =  ecoliinta.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/101
    eintadftoplot = eintadftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
eintadftoplot
eintadftoplot = eintadftoplot.round(2)
print(eintadftoplot)
Ecoliintatolineplot = eintadftoplot[['%ofpredictions','%ofTruePositives']]
Ecoliintatolineplot


import matplotlib.pyplot as plt
import numpy as np

plt.axis([10, 100, 10, 100])
plt.plot(Ecolirftolineplot[['%ofpredictions']],Ecolirftolineplot[['%ofTruePositives']],'g-',)
plt.show()


import matplotlib.pyplot as plt
import numpy as np

plt.axis([10, 100, 10, 100])
plt.plot(Ecoliintatolineplot[['%ofpredictions']],Ecoliintatolineplot[['%ofTruePositives']],'b-',)
plt.show()



import matplotlib.pyplot as plt
import numpy as np

plt.axis([10, 100, 10, 100])
plt.plot(Ecolicopratolineplot[['%ofpredictions']],Ecolicopratolineplot[['%ofTruePositives']],'r-')
plt.show()



import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

ax.set_xlim([10, 100])
ax.set_ylim([0, 110])

ax.set_xlabel('% of predictions')
ax.set_ylabel('% of true positives')

ax.plot(Ecolicopratolineplot[['%ofpredictions']],Ecolicopratolineplot[['%ofTruePositives']], 'r-', label='CopraRNA')
ax.plot(Ecolirftolineplot[['%ofpredictions']],Ecolirftolineplot[['%ofTruePositives']], 'g-', label='sRNARFTarget')
ax.plot(Ecoliintatolineplot[['%ofpredictions']],Ecoliintatolineplot[['%ofTruePositives']], 'b-', label='IntaRNA')


legend = ax.legend(loc='lower right', shadow=True)
ax.set_title('Ecoli')
plt.savefig('Ecoli_%_plot.pdf')
plt.show()



import pandas as pd
pccrf = pd.read_csv('./sRNARFTarget_PCC.txt', header = 0,sep = '\t')
pccrf.head()

pccrfrows = (pccrf.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
pccrfdftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((pccrfrows*i)/100)
    nows =  pccrf.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/20
    pccrfdftoplot = pccrfdftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
pccrfdftoplot
pccrfdftoplot = pccrfdftoplot.round(2)
print(pccrfdftoplot)
pccrftolineplot = pccrfdftoplot[['%ofpredictions','%ofTruePositives']]
pccrftolineplot




import pandas as pd
pcccopra = pd.read_csv('./CopraRNA_PCC.txt', header = 0,sep = '\t')
pcccopra.head()

pcccoprarows = (pcccopra.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
pcccopradftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((pcccoprarows*i)/100)
    nows =  pcccopra.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/20
    pcccopradftoplot = pcccopradftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
pcccopradftoplot
pcccopradftoplot = pcccopradftoplot.round(2)
print(pcccopradftoplot)
pcccopratolineplot = pcccopradftoplot[['%ofpredictions','%ofTruePositives']]
pcccopratolineplot




import pandas as pd
pccinta = pd.read_csv('./IntaRNA_PCC.txt', header = 0,sep = '\t')
pccinta.head()

pccintarows = (pccinta.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
pccintadftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((pccintarows*i)/100)
    nows =  pccinta.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/20
    pccintadftoplot = pccintadftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
pccintadftoplot
pccintadftoplot = pccintadftoplot.round(2)
print(pccintadftoplot)
pccintatolineplot = pccintadftoplot[['%ofpredictions','%ofTruePositives']]
pccintatolineplot



import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

ax.set_xlim([10, 100])
ax.set_ylim([0, 110])

ax.set_xlabel('% of predictions')
ax.set_ylabel('% of true positives')

ax.plot(pcccopratolineplot[['%ofpredictions']],pcccopratolineplot[['%ofTruePositives']], 'r-', label='CopraRNA')
ax.plot(pccrftolineplot[['%ofpredictions']],pccrftolineplot[['%ofTruePositives']], 'g-', label='sRNARFTarget')
ax.plot(pccintatolineplot[['%ofpredictions']],pccintatolineplot[['%ofTruePositives']], 'b-', label='IntaRNA')


legend = ax.legend(loc='lower right', shadow=True)
ax.set_title('PCC')
plt.savefig('PCC_%_plot.pdf')
plt.show()



import pandas as pd
multorf = pd.read_csv('./sRNARFTarget_Multocida.txt', header = 0,sep = '\t')
multorf.head()

multorfrows = (multorf.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
mrfdftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((multorfrows*i)/100)
    nows =  multorf.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/22
    mrfdftoplot = mrfdftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
mrfdftoplot
mrfdftoplot = mrfdftoplot.round(2)
print(mrfdftoplot)
multorftolineplot = mrfdftoplot[['%ofpredictions','%ofTruePositives']]
multorftolineplot




import pandas as pd
multocopra = pd.read_csv('./CopraRNA_Multocida.txt', header = 0,sep = '\t')
multocopra.head()

multocoprarows = (multocopra.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
mcopradftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((multocoprarows*i)/100)
    nows =  multocopra.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/22
    mcopradftoplot = mcopradftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
mcopradftoplot
mcopradftoplot = mcopradftoplot.round(2)
print(mcopradftoplot)
multocopratolineplot = mcopradftoplot[['%ofpredictions','%ofTruePositives']]
multocopratolineplot




import pandas as pd
multointa = pd.read_csv('./IntaRNA_Multocida.txt', header = 0,sep = '\t')
multointa.head()

multointarows = (multointa.shape)[0]

column_names = ["%ofpredictions", "No.ofpredictions", "No.ofTruePositives", "%ofTruePositives"]
mintadftoplot = pd.DataFrame(columns = column_names)

for i in range(10,110,10):
    perc = int((multointarows*i)/100)
    nows =  multointa.head(n=perc)
    TPdf = nows.loc[nows['Class'] == 1]
    TP = (TPdf.shape)[0]
    TPperc = (100*TP)/22
    mintadftoplot = mintadftoplot.append({'%ofpredictions' : int(i), 'No.ofpredictions' : perc , 'No.ofTruePositives' : TP,'%ofTruePositives' :TPperc } , ignore_index=True)
    
mintadftoplot
mintadftoplot = mintadftoplot.round(2)
print(mintadftoplot)
multointatolineplot = mintadftoplot[['%ofpredictions','%ofTruePositives']]
multointatolineplot




import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

ax.set_xlim([10, 100])
ax.set_ylim([0, 110])

ax.set_xlabel('% of predictions')
ax.set_ylabel('% of true positives')

ax.plot(multocopratolineplot[['%ofpredictions']],multocopratolineplot[['%ofTruePositives']], 'r-', label='CopraRNA')
ax.plot(multorftolineplot[['%ofpredictions']],multorftolineplot[['%ofTruePositives']], 'g-', label='sRNARFTarget')
ax.plot(multointatolineplot[['%ofpredictions']],multointatolineplot[['%ofTruePositives']], 'b-', label='IntaRNA')


legend = ax.legend(loc='lower right', shadow=True)
ax.set_title('Multocida')
plt.savefig('Multocida_%_plot.pdf')
plt.show()





