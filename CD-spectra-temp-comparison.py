#!/usr/bin/env ipython

# -*- coding: utf-8 -*-




import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib.ticker as ticker
import itertools
import collections 
import matplotlib as mpl
import matplotlib.artist as artists
from figure_dict import *

#get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')


# In[2]:


#input
fname = "input/CD-spectra/20170727-SsHelixMutant-Spectra.xlsx"
temp1 = "RT"
temp2 = "30"


#Output
dirOut = "output/"
experiment = "CD-spectra-temp-comparison"
fout = dirOut+experiment+'.pdf'
fout2 = dirOut+experiment+'-wtonly.pdf'
fout3 = dirOut+experiment+'-percentsignal.csv'


# In[3]:


xls = pd.ExcelFile(fname)
#mre for mutants: cols = wavelength, mut1, mut2...
df = xls.parse(temp1) 
df2 = xls.parse(temp2) 


x = df["Wavelength"]
x2 = df2["Wavelength"]
#remove wavelength, change to mdeg 
df =   df.iloc[:,1:]/1000
df2 = df2.iloc[:,1:]/1000
dict_temp1 = df.to_dict(orient = "list")
dict_temp2 = df2.to_dict(orient = "list")
protlist = list(df)
protlist2 = list(df2)


# In[4]:


#read in list of files, use to get proteins  
marker_p =   ["o"] * len(protlist)
markersize_p =  [40] * len(protlist)
markerdict = dict(zip(protlist, marker_p))
sizedict = dict(zip(protlist, markersize_p))

marker_p =   ["*"] * len(protlist2)
markersize_p =  [90] * len(protlist2)
markerdict2 = dict(zip(protlist2, marker_p))
sizedict2 = dict(zip(protlist2, markersize_p))



# In[5]:


### WT ONLY 
fig, ax = plt.subplots(figsize=(8,8))
 
fig_title = u'CD spectra - MRE vs wavelength' 
fig.suptitle(fig_title, fontsize=20)
plt.subplots_adjust(hspace=0.1)
fig_subtitle = "10mM KPi pH 7.2"

q = "SsWT"
ax.scatter(x, dict_temp1[q], c= colordict[q], marker=markerdict[q], s=sizedict[q],label=(q), linewidth='1', edgecolor='black')
ax.scatter(x, dict_temp2[q], c= colordict[q], marker=markerdict2[q], s=sizedict2[q],label=(q), linewidth='1', edgecolor='black')
#legend = ax.legend(loc="lower right", borderaxespad=1, ncol=2, handletextpad=0.1)        
str1 = u'\u26aa ' + temp1 + " "
str2 = u' \u2606 ' + str(temp2) + '$^\circ$C'
ax.text(254, -10.0, str1, ha='right', va='center')
ax.text(255.5, -10.5, str2, ha='right', va='center')
ax.set_xlim([215,260])
ax.set_ylim([-13,0])
ax.set_xlabel(labeldict["wavelength"])
ax.set_ylabel(labeldict["MRE-spectra"])

fig.tight_layout()
fig.subplots_adjust(top = 0.88)
#fig.show()
#fig.savefig(fout2)


# In[6]:


###--- plot scatter and line graph of urea titration for each mutant separately
fig, axarr = plt.subplots(3, 2, figsize=(15,15))

ax1 = axarr[0, 0]
ax2 =axarr[0, 1]
ax3 =axarr[1, 0]
ax4 =axarr[1, 1]
ax5 =axarr[2, 0]
ax6 =axarr[2, 1]

fig_title = u'CD spectra - MRE vs wavelength      ' 
fig.suptitle(fig_title, fontsize=20)
plt.subplots_adjust(hspace=0.1)
fig_subtitle = "10mM KPi pH 7.2"

counter = 0
for i in range(0,3):
    for j in range (0,2):
        counter = counter+1
        q = "SsWT"
        p = protlist[counter] ##mutant
        plabel = p.replace("_", "/")
        axarr[i, j].scatter(x, dict_temp1[p], c= colordict[p], marker=markerdict[p], s=sizedict[p], label=(plabel), linewidth='0.5', edgecolor='black')
        axarr[i, j].scatter(x, dict_temp1[q], c= colordict[q], marker=markerdict[q], s=sizedict[q],label=(q), linewidth='1', edgecolor='black')
        axarr[i, j].scatter(x, dict_temp2[p], c= colordict[p], marker=markerdict2[p], s=sizedict2[p], label=(plabel), linewidth='0.5', edgecolor='black')
        axarr[i, j].scatter(x, dict_temp2[q], c= colordict[q], marker=markerdict2[q], s=sizedict2[q],label=(q), linewidth='1', edgecolor='black')
        legend = axarr[i, j].legend(loc="lower right", borderaxespad=1, ncol=2, handletextpad=0.1)        
        str1 = u'\u26aa ' + temp1 + " "
        str2 = u' \u2606 ' + str(temp2) + '$^\circ$C'
        if (p == "I45A_S70A" or p == "I45A_M73A"): 
            axarr[i, j].text(236, -9.3, str1, ha='center', va='center')
            axarr[i, j].text(250, -9.3, str2, ha='center', va='center')
        else:
            axarr[i, j].text(244, -9.3, str1, ha='center', va='center')
            axarr[i, j].text(254, -9.3, str2, ha='center', va='center')
        axarr[i, j].set_xlim([215,260])
        axarr[i, j].set_ylim([-13,0])
        axarr[i, j].set_xlabel(u'Wavelength (nm)')
        axarr[i, j].set_ylabel(r'$MRE (10^3 deg \cdot cm^2 \cdot dmol^{-1})$')

fig.tight_layout()
fig.subplots_adjust(top = 0.88)
#fig.show()
fig.savefig(fout2)


# In[7]:


xls = pd.ExcelFile(fname)
#mre for mutants: cols = wavelength, mut1, mut2...
df = xls.parse(temp1) 
df2 = xls.parse(temp2) 


# In[8]:


RT = df.loc[df['Wavelength'] == 222]
thirty = df2.loc[df['Wavelength'] == 222]
percentsignal = thirty.div(RT, axis = "row")*100
#df.at['Wavelength', 'Wavelength'] = 10
percentsignal.set_index('Wavelength')
percentsignal.index = ["222"]
percentsignal = percentsignal.drop(['Wavelength'], axis=1).round(decimals=2)
percentsignal.to_csv(fout3, index=True, index_label = "Wavelength")
#percentsignal

