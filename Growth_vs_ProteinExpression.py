
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
##print(plt.style.available)


# In[2]:


dirIn = "input/GrowthAssay/"
dirOut = "output/"
fname = "SsMutant"
extracaption = "-unlabeled"
temp = 30 
fin =  dirIn+fname+".xlsx"
fout = fname+extracaption+'-selection_coef_vs_protein_conc.pdf'

xls = pd.ExcelFile(fin)
#cols = mutant, selection coef
df = xls.parse("western") 


# In[3]:


df["color_p"] = df["Protein"].map(colordict)
df["label_p"] = df["Protein"].str.replace('_','\n')

#drop vector and any other rows with missing values 
df = df.dropna()
df


# In[4]:


def autolabel(p, x, y, ax):
    for prot, x_val, y_val in zip(p, x, y):
        invertlabel = ["SsWT", "I107A"]
        moveleftlabel = ["I45A\nI107A"]
        moverightlabel = ["I45A\nS70A"]
        if (prot in invertlabel): 
            ax.text(x_val, y_val-0.02, prot, ha='center', va='top')
        elif (prot in moveleftlabel): 
            ax.text(x_val-0.1, y_val+0.01, prot, ha='center', va='bottom')
        elif (prot in moverightlabel): 
            ax.text(x_val+0.15, y_val-0.05, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val+0.01, prot, ha='center', va='bottom')
    return

def scatter(ax, x, y, col_p, p, x_col, y_col):
    wt_x = df.loc[df["Protein"] == "SsWT", x_col].item()
    wt_y = df.loc[df["Protein"] == "SsWT", y_col].item()
    ax.axvline(x=wt_x, color = "grey", linewidth='1', linestyle = ":", zorder=0)
    ax.axhline(y=wt_y, color = "grey", linewidth='1', linestyle = ":", zorder=0)
    ax.scatter(x, y, c=col_p , marker="o", s=50, linewidth='0.5', edgecolor='black', label='_nolegend_')
    #autolabel(p, x, y, ax)
    padding = 0.1
    ax.set_xlim([min(x)-padding,max(x)+padding])
    ax.set_ylim([min(y)-padding,max(y)+padding])
    #plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    ax.tick_params(axis='both') 
    return


# In[5]:


x_col = "ProteinExpression"
y_col = "Selcoeff_all"

x_label = labeldict[x_col]
y_label = labeldict[y_col]

#x_ax = axesdict[x_col]
#y_ax = axesdict[y_col]

fig = plt.figure(figsize=(8,8)) 
ax1 = plt.subplot()
fig_title = y_label + " vs " + x_label
fig.suptitle(fig_title)
plt.subplots_adjust(hspace=0.1)
ax1.set_xlabel(x_label)
ax1.set_ylabel(y_label)

x = df[x_col]
y = df[y_col]
col_p = df["color_p"]
p = df["label_p"]

scatter(ax1, x, y, col_p, p, x_col, y_col)


fig.tight_layout()
fig.subplots_adjust(top = 0.88)
fig.show()


# In[6]:


fig.savefig(dirOut+fout)

