
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
fout = dirOut+"SsMutant-30C-ddG-activation-energy.pdf"
df = pd.read_csv(dirOut+fname+'.csv', index_col = "Protein")
df = df.drop('SsWT', axis=0)
temp = 30
extracaption = "-"+str(temp)+"C"
df['Prot'] = df.index
df


# In[3]:


def autolabel(rects, ddg, labels):
    for rect, h in zip(rects, ddg):
        height = h
        if (height >= 0):
            ax.text(rect.get_x() + rect.get_width()/2, height+.03, '%.2f' % height, ha='center', va='center')
        if (height < 0):
            ax.text(rect.get_x() + rect.get_width()/2, height-.03, '%.2f' % height, ha='center', va='center')


# In[4]:


#ddg-kcat (kcal/mol)	ddg-Keff (kcal/mol)

#### 
ddg = df[["ddg-kcat (kcal/mol)", "ddg-Keff (kcal/mol)"]]
labels = df.index.str.replace('_' , '\n')

fig = plt.figure(figsize=(12, 8)) 
gs = gridspec.GridSpec(nrows=1, ncols=1)
ax = plt.subplot(gs[0])


N = df.shape[0]
ind = np.arange(N)  # the x locations for the groups
width = 0.4 # the width of the bars

rects3 = ax.bar(ind, ddg["ddg-kcat (kcal/mol)"], width, color="white",edgecolor='black')
rects4 = ax.bar(ind+width, ddg["ddg-Keff (kcal/mol)"] , width, color="white", edgecolor='black', hatch='//')

rects1 = ax.bar(ind, ddg["ddg-kcat (kcal/mol)"] , width, color=[colordict[p] for p in df["Prot"]], edgecolor='black')
rects2 = ax.bar(ind+width, ddg["ddg-Keff (kcal/mol)"] , width, color=[colordict[p] for p in df["Prot"]], edgecolor='black', hatch='//')
#df["ddg-EnergyActivation (kcal/mol)"].plot(kind='bar', color=[colordict[p] for p in df["Prot"]])

fig_title = u'$\u0394\u0394G_{T}^{\u2021}$'+ ' transition state stabilization free energy'
sub_title = u'$\u0394\u0394G_{T}^{\u2021}$' + u'$=-RTln[(kcat/Km)_{Mut}/(kcat/Km)_{WT}]$'
plt.suptitle(fig_title)
#ax.set_title(sub_title, fontsize=12)
ax.set_ylabel(u'$\u0394\u0394G^{\u2021}$ ' + r'$(cal \cdot mol^{-1})$')
ax.set_xticks(ind + 2*width / 2)
ax.set_xticklabels(labels)
ax.set_ylim([min(ddg["ddg-Keff (kcal/mol)"])*1.2, max(ddg["ddg-Keff (kcal/mol)"])*1.2])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.legend((rects3[0], rects4[0]), (u'$\u0394\u0394G_{Kcat}$', u'$\u0394\u0394G_{Keff}$'), loc=2)
ax.margins(0.03)
autolabel(rects1, ddg["ddg-kcat (kcal/mol)"], labels)
autolabel(rects2, ddg["ddg-Keff (kcal/mol)"], labels)
ax.tick_params(axis='both') 
ax.axhline(y=0, color = "black")
fig.show()


# In[5]:


fig.savefig(fout)

