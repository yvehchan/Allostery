
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


fout = 'output/ddGNI_vs_ddGkcat.pdf'

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
df = df.drop("SsWT", axis=0)

a = [["ddG_NI", "ddG_total"],["ddg-kcat (kcal/mol)","ddg-Keff (kcal/mol)"]]
xyzcombo = list(itertools.product(*a))
df[["ddG_NI", "ddG_total","ddg-kcat (kcal/mol)","ddg-Keff (kcal/mol)"]]


# In[3]:


def ceil_power_of_10(n):
    exp = math.log(abs(n), 10)
    exp = math.ceil(exp)
    return 10**exp

def autolabel(p, x, y, ax, intercept, slope, minx, maxx):
    for prot, x_val, y_val in zip(p, x, y):
        ypred = intercept+slope*x_val
        y_shift = ceil_power_of_10(min(y)) * 0.02
        x_shift = 0.15
        if x_val+x_shift > maxx:
            x_val = x_val - x_shift
        if x_val-x_shift < minx:
            x_val = x_val + x_shift
        ha = "center"
    
        if (y_val < ypred or prot == "I45K" or prot == "I45K"):
            va = "top"
            if prot == "I45K":
                ha = "right"
                va = "bottom"
            if prot == "I45A":
                ha = "left"
                va = "top"
            ax.text(x_val, y_val-y_shift, prot, ha=ha, va=va, fontsize=13)
        else:
            ax.text(x_val, y_val+y_shift, prot, ha=ha, va="bottom", fontsize=13)
    return


# In[4]:


def get2DScatter (df, plotnum, cols, labels ): 
    x = df[cols[0]]
    y = df[cols[1]]
    ax = fig.add_subplot(2, 2, plotnum)
    
    #x = x.sort_values() 
    X = np.vander(x, 2)
    model1 = sm.OLS(y,X)
    results1 = sm.OLS(y,X).fit()
    xx = np.linspace(min(x), max(x), 100)
    XX = np.vander(xx, 2)
    ax.plot(xx, results1.predict(XX), 'k--', label='Poly n=1, $R$=%.2f' % math.sqrt( results1.rsquared), alpha=0.9, zorder = 0)

    # 2-nd order polynomial
    # add a constant and quadratic term
    X2 = np.vander(x, 3)
    model2 = sm.OLS(y,X2)
    results2 = sm.OLS(y,X2).fit()
    XX = np.vander(xx, 3)
    ax.plot(xx, results2.predict(XX), 'g--', label='Poly n=2, $R$=%.2f' % math.sqrt( results2.rsquared), alpha=0.9, zorder = 0)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    #wt_x = df.get_value("SsWT", cols[0])
    #wt_y = df.get_value("SsWT", cols[1])
    #ax.axvline(x=wt_x, color = "lightgrey", linestyle = "dashed")
    #ax.axhline(y=wt_y, color = "lightgrey", linestyle = "dashed")
    ax.plot(x, intercept+slope*x, '--k', label='_nolegend_')
    ax.scatter(x, y, c= df["color_p"] , marker="o", s=30, linewidth='0.5', edgecolor='black', label='_nolegend_')
    autolabel(df["label_p"], x, y, ax, intercept, slope, min(x), max(x))
    #if (slope > 0):
    #    ax.text(max(x), min(y), u'$R= %.2f$' % r_value, ha='left', va='center', fontsize=16, color = "black")
    #else:
    #ax.text(max(x), min(y), u'$R= %.2f$' % r_value, ha='right', va='center', fontsize=16, color = "black")	
    #x_shift = ceil_power_of_10(min(x)) * 0.1
    #y_shift = ceil_power_of_10(min(y)) * 0.2
    x_shift = 0.3
    y_shift = 0.3
    ax.set_xlim([min(x)-x_shift,max(x)+x_shift])	
    ax.set_ylim([min(y)-y_shift,max(y)+2*y_shift])
    fig_subtitle = labels[0] + " " + labels[1]
    #ax.set_title(fig_subtitle, fontsize=8)
    ax.set_xlabel(labels[0], fontsize = 18)
    ax.set_ylabel(labels[1], fontsize = 20)
    ax.tick_params(axis='both', labelsize=16) 
    ax.legend(loc = "best", fontsize=16)

    return 

fig = plt.figure(figsize=(12, 12))

for i, xyz in enumerate(xyzcombo):
	x_col = str(xyz[0])
	y_col = str(xyz[1])
	x_label = labeldict[x_col]
	y_label = labeldict[y_col]

	cols = [x_col, y_col]
	labels = [x_label, y_label]

	get2DScatter (df, i+1, cols, labels)

fig_title = u'ddG thermodynamic stability vs ddG transition state' 
fig.suptitle(fig_title, fontsize=12)
fig.tight_layout()	
#plt.subplots_adjust(hspace=0.1)	
fig.subplots_adjust(top = 0.90, hspace=0.2)
fig.show()


# In[5]:


fig.savefig(fout)

