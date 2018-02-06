#!/usr/bin/env ipython
# -*- coding: utf-8 -*-




import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from figure_dict import *

#get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib





#input file - thermal melt first derivative	
fin = "input/GrowthAssay/SsMutant.xlsx"
xls = pd.ExcelFile(fin)
#cols = mutant, selection coef
df = xls.parse("singles_clamps") 

#Output
dirOut = "output/"
experiment = "fitness"
fout = dirOut+experiment+'-avg-scatter.pdf'
csvout = dirOut+experiment+'-avg-scatter.csv'

colnames = list(df)
dfmax = df.max().max()
dfmin = df.min().min()
df["rep_num"] = ["r1", "r2", "r3"]

#print (fout)
#df


# In[3]:


df.describe().to_csv(csvout, index=True)


# In[4]:


def getStripPlot (double, plotnum): 
    ax = fig.add_subplot(2, 2, plotnum)
    mut_1 = double.split("_")[0]
    mut_2 = double.split("_")[1]
    cur_df = df.filter(items=[mut_1, mut_2, double, "rep_num"])
    df2 = cur_df.melt('rep_num').dropna()
    #cur_df.mean().round(2)
    cur_df.rename(columns = {double:mut_1+"\n"+mut_2}, inplace = True)
    meanpointprops = dict(marker='*', linewidth=3, markeredgecolor='black', markerfacecolor="None")
    medianprops = dict(linewidth=0)
    meanlineprops = dict(linestyle='--', linewidth=2.5, color='black')
    colordict = {"r1": "lightgray", "r2": "darkgray", "r3": "black"}
    sns.swarmplot(x="variable", y="value", hue="rep_num", palette=colordict, edgecolor="k",linewidth=1, dodge=True, size=10, data=df2)    
    sns.boxplot(data=cur_df, width = 0.4,
                meanprops=meanlineprops, meanline=True, showmeans=True,
                medianprops=medianprops, showcaps=False,showbox=False,
                boxprops={'facecolor':'None'},showfliers=False,whiskerprops={'linewidth':0})
    ax.set_ylabel(u's') 
    ax.set_xlabel(u'Mutant')
    ax.set_title(mut_1+" "+mut_2)
    ax.set_ylim([dfmin-0.1, dfmax+0.1])
    ax.tick_params(axis='x',which='both',bottom='off',top='off')
    if (plotnum == 1):
        ax.legend(bbox_to_anchor=(0.02, 0.03), loc="lower left", borderaxespad=0.)
    else: 
        ax.legend_.remove()
    ax.axhline(y=0, color = "black")
    return


# In[5]:


fig = plt.figure(figsize=(12, 8))
doubles = ["I45A_S70A","I45A_M73A","I45A_I107A","I107A_D128A"]
for index, double in enumerate(doubles):
    getStripPlot(double, index+1)
fig.tight_layout()	
#fig.show()


# In[6]:


fig.savefig(fout)

