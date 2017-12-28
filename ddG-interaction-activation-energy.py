
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
from matplotlib.offsetbox import AnchoredText
import matplotlib.artist as artists
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from figure_dict import *
import statsmodels.api as sm
import statsmodels.formula.api as smf
import lmfit as lm
import re
from scipy import stats
from scipy.optimize import curve_fit
import itertools
from collections import defaultdict
import numbers
from uncertainties import ufloat
from uncertainties import umath
from uncertainties.umath import *
from uncertainties import unumpy
from decimal import Decimal
import math

get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib
##print(plt.style.available)


# In[2]:


dirIn = 'input/'
dirOut = 'output/'

fin = "output/SsMutant-double-pooled-30C-initialvelocity-mm.csv"
fout = dirOut+"ddG-activation-energy.pdf"
df = pd.read_csv(fin, index_col = "Protein")
df['Prot'] = df.index
#df[df['Prot'].isin(["I45A_M73A", "I45A_S70A"])]
#print (df)
int_kcat_I45A_S70A = df["ddg-kcat (kcal/mol)"]["I45A_S70A"] - (df["ddg-kcat (kcal/mol)"]["I45A"] + df["ddg-kcat (kcal/mol)"]["S70A"])
int_kcat_I45A_M73A = df["ddg-kcat (kcal/mol)"]["I45A_M73A"] - (df["ddg-kcat (kcal/mol)"]["I45A"] + df["ddg-kcat (kcal/mol)"]["M73A"])
int_keff_I45A_S70A = df["ddg-Keff (kcal/mol)"]["I45A_S70A"] - (df["ddg-Keff (kcal/mol)"]["I45A"] + df["ddg-Keff (kcal/mol)"]["S70A"])
int_keff_I45A_M73A = df["ddg-Keff (kcal/mol)"]["I45A_M73A"] - (df["ddg-Keff (kcal/mol)"]["I45A"] + df["ddg-Keff (kcal/mol)"]["M73A"])

ddgni_list_values =  [int_kcat_I45A_S70A, int_kcat_I45A_M73A]
#ddgni_err_list_values =  add['interaction_NI_err'].tolist()
ddgiu_list_values =  [int_keff_I45A_S70A, int_keff_I45A_M73A]
#ddgiu_err_list_values =  add['interaction_IU_err'].tolist()

df


# In[3]:


minvalue = min(min(ddgni_list_values), min(ddgiu_list_values))
maxvalue = max(max(ddgni_list_values), max(ddgiu_list_values))
protlist = ["I45A_S70A", "I45A_M73A"]
labels = [w.replace('_', '\n') for w in protlist]

def autolabel(rects, list_values):
	for rect, h in zip(rects, list_values):
		height = h
		if (height >= 0):
			ax.text(rect.get_x() + rect.get_width()/2, height+0.004, '%.2f' % height, 
                    ha='center', va='center')
		if (height < 0):
			ax.text(rect.get_x() + rect.get_width()/2, height-0.004, '%.2f' % height, 
                    ha='center', va='center')


# In[4]:


colorlist = ["green", "blue"]
#mpl.rcParams['hatch.linewidth'] = 3

N = 2
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
fig = plt.figure(figsize=(12, 8)) 
plt.rcParams["patch.force_edgecolor"] = True
gs = gridspec.GridSpec(nrows=1, ncols=1)
ax = plt.subplot(gs[0])

##dummy plot for legend marker
rects3 = ax.bar(ind, ddgni_list_values, width, color="white",edgecolor='black')
rects4 = ax.bar(ind + width, ddgiu_list_values, width, color="white", edgecolor='black', hatch='//')


rects1 = ax.bar(ind, ddgni_list_values, width, color=colorlist,edgecolor='black')
rects2 = ax.bar(ind+width, ddgiu_list_values, width, color=colorlist)
rects5 = ax.bar(ind+width, ddgiu_list_values, width, color="None", hatch='//')

# add some text for labels, title and axes ticks
ax.set_ylabel(u'\u03B4 ' + r'$(cal \cdot mol^{-1})$')
ax.set_title(u'\u03B4 ' + r'$(cal \cdot mol^{-1})$' + ' interaction energy of double mutant')
ax.set_xticks(ind + width )
ax.tick_params(axis='both') 
ax.set_xticklabels(labels)

ax.set_ylim([minvalue*1.1, 0])
autolabel(rects1, ddgni_list_values)
autolabel(rects2, ddgiu_list_values)

ax.legend((rects3[0], rects4[0]), (u'$\u03B4_{kcat}$', u'$\u03B4_{keff}$'),  prop={'size': 20}, loc=4)
ax.axhline(y=0, color = "black")
ax.margins(0.03)


# In[5]:


fig.savefig(fout)

