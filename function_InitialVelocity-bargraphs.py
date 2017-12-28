
# coding: utf-8

# In[1]:


import csv
import pandas as pd
import numpy as np
from numpy.random import normal

##matplotlib.use('TkAgg')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
import seaborn as sns
from figure_dict import *
import statsmodels.api as sm
import lmfit as lm
import re
from scipy import stats
from scipy.optimize import curve_fit
from collections import defaultdict
import numbers
from uncertainties import ufloat
from uncertainties import umath
from uncertainties.umath import *
from uncertainties import unumpy
from decimal import Decimal

get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib
##print(plt.style.available)


# In[2]:


dirOut = 'output/'
fname = "SsMutant-double-pooled-30C-initialvelocity-mm"
fin = dirOut+fname+".csv"
fout = dirOut+"SsMutant-30C-initialvelocity-mm-bargraph.pdf"
cols = ["Protein","conc", "vmax", "vmaxerr", "kcat", "kcaterr", "km", "kmerr", "kcatkm", "kcatkmerr", "kcatM","ddg-kcat", "ddg-Keff", "kcal", "Keff"]
df = pd.read_csv(fin, names=cols, skiprows =1)
df["Prot"] = df["Protein"]
df.set_index('Protein', inplace=True)
df


# In[3]:


orderdict = {
'SsWT': 1, 
'I45A': 2, 
'I45K': 3, 
'S70A': 4,
'M73A': 5, 
'M73I': 6, 
'I45A_S70A': 7, 
'I45A_M73A': 8,
'I107A': 9, 
'I107K': 10,
'I45A_I107A': 11,
'D128A': 12, 
'I107A_D128A': 13, 
'D121A': 14, 
'K207A': 15, 
'D121A_K207A': 16, 
'D165A': 17
}

cols = ["Protein", "conc", "vmax", "vmaxerr", "kcat", "kcaterr", "km", "kmerr", "kcatkm", "kcatkmerr", "ddg-kcat", "ddg-Keff", "kcal", "Keff"]

def autolabel(rects, ddg, ax, labels, scale, numsig):
    for rect, h in zip(rects, ddg):
        height = h
        if numsig == 2:
            s = '%.1f'
        elif numsig == 0:
            s = '%d'
        else:
            s = '%.3f'

        if (height >= 0):
            ax.text(rect.get_x() + rect.get_width()/2, height+scale, s % height, 
                    ha='center', va='center', fontsize=16)
        if (height < 0):
            ax.text(rect.get_x() + rect.get_width()/2, height-scale, s % height, 
                    ha='center', va='center', fontsize=16)
            #ax.text(rect.get_x() + rect.get_width()/2, 0.02, p , ha='center', va='bottom', color = p_color)


# In[4]:


vmax = df["vmax"]
vmaxerr = df["vmaxerr"]
kcat = df["kcat"]
kcaterr =df["kcaterr"]
km = df["km"]
kmerr = df["kmerr"]
kcatkm = df["kcatkm"]
kcatkmerr = df["kcatkmerr"]
labels = df.index.str.replace('_' , '\n')


# In[5]:


fig = plt.figure(figsize=(12, 12)) 
gs = gridspec.GridSpec(nrows=3, ncols=1)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])

N = df.shape[0]
ind = np.arange(N)  # the x locations for the groups
width = 0.4       # the width of the bars

ax0.axhline(y=kcat["SsWT"], color = "grey", linestyle = ":", zorder=0)
ax1.axhline(y=km["SsWT"], color = "grey", linestyle = ":", zorder=0)
ax2.axhline(y=kcatkm["SsWT"], color = "grey", linestyle = ":", zorder=0)


rects0 = ax0.bar(ind, kcat, width, yerr=kcaterr, color=[colordict[p] for p in df["Prot"]], edgecolor='black')
rects1 = ax1.bar(ind, km ,  width, yerr=kmerr,color=[colordict[p] for p in df["Prot"]], edgecolor='black')
rects2 = ax2.bar(ind, kcatkm, width, yerr=kcatkmerr, color=[colordict[p] for p in df["Prot"]], edgecolor='black')

fig_title = u'Kinetic parameters'
#sub_title = u'$\u0394\u0394G_{T}^{\u2021}$' + u'$=-RTln[(kcat/Km)_{Mut}/(kcat/Km)_{WT}]$'
plt.suptitle(fig_title, fontsize=14)
#ax.set_title(sub_title, fontsize=12)
ax0.set_ylabel(u'$k_{cat}$' + u'($s^{-1})$')
ax1.set_ylabel(u'$K_{m}$' + u'(\u00B5M)')
ax2.set_ylabel(u'$k_{eff}$'+ u'(\u00B5M' + r'$\cdot s^{-1})$')

ax0.set_xticks(ind+width/2)
ax1.set_xticks(ind+width/2)
ax2.set_xticks(ind+width/2)

ax0.set_xticklabels(labels)
ax1.set_xticklabels(labels)
ax2.set_xticklabels(labels)

ax0.set_ylim([0, max(kcat)*1.2])
ax1.set_ylim([0, max(km)*1.2])
ax2.set_ylim([0, max(kcatkm)*1.2])

ax0.tick_params(axis='both') 
ax1.tick_params(axis='both') 
ax2.tick_params(axis='both') 

ax0.margins(0.03)
ax1.margins(0.03)
ax2.margins(0.03)

ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 

autolabel(rects0, kcat, ax0,  labels, 1.3, 2)
autolabel(rects1, km, ax1, labels, 60, 0)
autolabel(rects2, kcatkm, ax2, labels, 0.002, 3)

#ax.axhline(y=0, color = "black")
fig.show()


# In[6]:


fig.savefig(fout)

