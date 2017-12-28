
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


fout = 'output/dG_vs_kcat.pdf'

fin = "output/SsMutant-ddG.csv"
dG = pd.read_csv(fin)

fin = "input/GrowthAssay/SsMutant.xlsx"
xls = pd.ExcelFile(fin)
dbl = xls.parse("selcoeff")

fin = "output/SsMutant-double-pooled-30C-initialvelocity-mm.csv"
kcat_temp = "30"
kcat = pd.read_csv(fin)


kcat["color_p"] = kcat["Protein"].map(colordict)
kcat["label_p"] = kcat["Protein"].str.replace('_','\n')


df = dbl.set_index('Protein').join(dG.set_index('Protein'), how = "inner")
df = df.join(kcat.set_index('Protein'), how = "inner")

df



# In[3]:


x_col = "dG_NI"
x2_col = "dG_total"
y_col = "kcat (1/s)"

x_label = labeldict[x_col]
x2_label = labeldict[x2_col]
y_label = labeldict[y_col]

x_ax = axesdict[x_col]
x2_ax = axesdict[x2_col]
y_ax = axesdict[y_col]
    


# In[4]:


def autolabel(p, x, y, ax, intercept, slope, minx, maxx):
    for prot, x_val, y_val in zip(p, x, y):
        ypred = intercept+slope*x_val
        x_shift = 0.1
        if x_val+x_shift > maxx:
            x_val = x_val - x_shift
        if x_val-x_shift < minx:
            x_val = x_val + x_shift
        y_shift = 0.3
        if (y_val > ypred):
            ax.text(x_val, y_val+y_shift, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val-y_shift, prot, ha='center', va='top')
    return


# In[5]:


def scatter(df, ax, x_col, y_col):
    wt_x = df.get_value("SsWT", x_col)
    wt_y = df.get_value("SsWT", y_col)
    
    temp = df.sort_values(by = [x_col])
    x = temp[x_col]
    y = temp[y_col]
    col_p = temp["color_p"]
    lab_p = temp["label_p"]

    ax.axvline(x=wt_x, color = "grey", linestyle = ":", zorder=0)
    ax.axhline(y=wt_y, color = "grey", linestyle = ":", zorder=0)
  
    ax.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black', label='_nolegend_')
    #ax.plot(x, intercept+slope*x, '--k')
    
    # 1-st order polynomial
    # add a constant row for the intercept
    #x = x.sort_values() 
    X = np.vander(x, 2)
    model1 = sm.OLS(y,X)
    results1 = sm.OLS(y,X).fit()
    xx = np.linspace(min(x), max(x), 100)
    XX = np.vander(xx, 2)
    ax.plot(xx, results1.predict(XX), 'k--', label='Poly n=1, $R$=%.2f' % math.sqrt( results1.rsquared), alpha=0.9)

    # 2-nd order polynomial
    # add a constant and quadratic term
    X2 = np.vander(x, 3)
    model2 = sm.OLS(y,X2)
    results2 = sm.OLS(y,X2).fit()
    XX = np.vander(xx, 3)
    ax.plot(xx, results2.predict(XX), 'g--', label='Poly n=2, $R$=%.2f' % math.sqrt( results2.rsquared), alpha=0.9)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    autolabel(lab_p, x, y, ax, intercept, slope, min(x), max(x))
    paddingx = 0.3
    paddingy = 1.5

    ax.set_xlim([min(x)-paddingx,max(x)+paddingx])
    ax.set_ylim([min(y)-paddingy,max(y)+paddingy])

    ax.legend(loc = 1, fontsize=16)
    return


# In[6]:


fig = plt.figure(figsize=(16, 8)) 
gs = gridspec.GridSpec(nrows=1, ncols=2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#fig = plt.figure(figsize=(8,8)) 
#ax1 = plt.subplot()
fig_title = u'dG vs kcat       ' 
fig.suptitle(fig_title, fontsize=20)
plt.subplots_adjust(hspace=0.1)
fig_subtitle = y_label + " vs " + x_label
ax1.set_title(fig_subtitle, fontsize=18)
ax1.set_xlabel(x_label, fontsize=18)
ax1.set_ylabel(y_label, fontsize=16)
fig_subtitle = y_label + " vs " + x2_label
ax2.set_title(fig_subtitle, fontsize=18)
ax2.set_xlabel(x2_label, fontsize=18)
ax2.set_ylabel(y_label, fontsize=16)

df1 = df.sort_values(x_col)
x = df1[x_col]
y = df1[y_col]
col_p = df1["color_p"]
p = df1["label_p"]

#print df[[x2_col,y_col,"color_p", "label_p"]]
df2 = df.sort_values(x2_col)
#print df2[[x2_col,y_col,"color_p","label_p"]]
x2 = df2[x2_col]
y2 = df2[y_col]
col_p2 = df2["color_p"]
p2 = df2["label_p"]

scatter(df, ax1, x_col, y_col)
scatter(df, ax2, x2_col, y_col)

fig.tight_layout()
fig.subplots_adjust(top = 0.88)
fig.show()


# In[7]:


fig.savefig(fout)

