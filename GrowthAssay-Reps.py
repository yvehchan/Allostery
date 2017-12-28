
# coding: utf-8

# In[1]:


import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from figure_dict import *

get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib


# In[2]:


dirIn = "input/GrowthAssay/"
dirOut = "output/"
fname = "SsMutant"
extracaption = ""
temp = 30 
fin =  dirIn+fname+".xlsx"
fout = fname+extracaption+'-selection_coef_boxplot.pdf'

xls = pd.ExcelFile(fin)
#cols = mutant, selection coef
df = xls.parse("singles_clamps") 
#df = df.pivot(index='Date', columns='Mutation', values='S')
colnames = list(df)
dfmax = df.max().max()
dfmin = df.min().min()


# In[3]:


def getStripPlot (double, plotnum): 
    ax = fig.add_subplot(2, 2, plotnum)
    mut_1 = double.split("_")[0]
    mut_2 = double.split("_")[1]
    print mut_1, mut_2, double
    cur_df = df.filter(items=[mut_1, mut_2, double])

    colorlist = []
    [colorlist.append(colordict[x]) for x in list(cur_df)]
    cur_df.rename(columns = {double:mut_1+"\n"+mut_2}, inplace = True)
    meanpointprops = dict(marker='*', linewidth=3, markeredgecolor='black', markerfacecolor="None", markersize=20)
    medianprops = dict(linewidth=0)
    meanlineprops = dict(linestyle='--', linewidth=2.5, color='black')
    sns.swarmplot(data=cur_df, size=12, palette=colorlist, edgecolor="k",linewidth=2)
    sns.boxplot(data=cur_df, width = 0.4,
                meanprops=meanlineprops, meanline=True, showmeans=True,
                medianprops=medianprops, showcaps=False,showbox=False,
                boxprops={'facecolor':'None'},showfliers=False,whiskerprops={'linewidth':0})

    ax.set_ylabel(u's') 
    ax.set_xlabel(u'Mutant')
    ax.set_title(mut_1+" "+mut_2)
    ax.set_ylim([dfmin-0.1, dfmax+0.1])
    ax.tick_params(axis='x',which='both',bottom='off',top='off')
    #ax.set_xticklabels([])
    #ax.set_xticks([])

    ax.axhline(y=0, color = "black")
    return

fig = plt.figure(figsize=(12, 8))
doubles = ["I45A_S70A","I45A_M73A","I45A_I107A","I107A_D128A"]
for index, double in enumerate(doubles):
    getStripPlot(double, index+1)
fig.tight_layout()	
fig.show()





# In[4]:


fig.savefig(dirOut+fout)

