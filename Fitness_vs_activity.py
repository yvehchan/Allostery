
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


def ceil_power_of_10(n):
    exp = math.log(abs(n), 10)
    exp = math.ceil(exp)
    return 10**exp

def autolabel(p, x, y, ax, intercept, slope):
	for prot, x_val, y_val in zip(p, x, y):
		ypred = intercept+slope*x_val
		y_shift = ceil_power_of_10(min(y)) * 0.01
		if (y_val >= ypred):
			ax.text(x_val, y_val+y_shift, prot, ha='center', va='bottom')
		else: 
			ax.text(x_val, y_val-y_shift, prot, ha='center', va='top')
	return
	
def getScatter (df, ax, x_col, y_col ): 
    temp = df.sort_values(by = [x_col])
    x = temp[x_col]
    y = temp[y_col]
    col_p = temp["color_p"]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    ax.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black')
    ax.plot(x, intercept+slope*x, '--k')
    #autolabel(df["label_p"], x, y, ax, intercept, slope)
    if (slope > 0):
        ax.text(min(x), max(y), u'$R= %.2f$' % r_value, ha='left', va='center', color = "black", fontsize=16)
    else:
        ax.text(max(x), max(y), u'$R= %.2f$' % r_value, ha='right', va='center', color = "black", fontsize=16)			
    x_shift = ceil_power_of_10(min(x)) * 0.1
    y_shift = ceil_power_of_10(min(y)) * 0.1
    ax.set_xlim([min(x)-x_shift,max(x)+x_shift])	
    ax.set_ylim([min(y)-y_shift,max(y)+y_shift])
    ax.axvline(x=x["SsWT"], color = "grey", linestyle = ":", zorder=0)
    ax.axhline(y=y["SsWT"], color = "grey", linestyle = ":", zorder=0)
    return


# In[3]:


fin = "input/GrowthAssay/SsMutant.xlsx"
fout = 'fitness_vs_activity'
selcoef_temp = "30"
xls = pd.ExcelFile(fin)
dbl = xls.parse("selcoeff")
 
fin = "output/SsMutant-double-pooled-30C-initialvelocity-mm.csv"
kcat_temp = "30"
kcat = pd.read_csv(fin)
kcat["color_p"] = kcat["Protein"].map(colordict)
kcat["label_p"] = kcat["Protein"].str.replace('_','\n')
df = dbl.set_index('Protein').join(kcat.set_index('Protein'), how = "inner")
df


# In[4]:



fig = plt.figure(figsize=(12, 12)) 
gs = gridspec.GridSpec(nrows=2, ncols=2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#ax3 = plt.subplot(gs[2])
#ax4 = plt.subplot(gs[3])

fig_title = u'Fitness vs enzyme activity' 
fig.suptitle(fig_title)
plt.subplots_adjust(hspace=0.1)

fig_subtitle = 'selection coefficient vs kcat'
ax1.set_title(fig_subtitle)
ax1.set_xlabel(u'$k_{cat}$'+ u' $(s^{-1})$')
ax1.set_ylabel('s')
fig_subtitle = 'Selection coefficient vs keff'
ax2.set_title(fig_subtitle)
ax2.set_xlabel(u'K$_{eff}$'+ u' $(\u00B5M^{-1} \cdot s^{-1})$')
ax2.set_ylabel('s')
#fig_subtitle = 'doubling time vs kcat'
#ax3.set_title(fig_subtitle)
#ax3.set_xlabel(u'$k_{cat} (s^{-1})$')
#ax3.set_ylabel('Doubling time (Hrs)')
#fig_subtitle = 'doubling time vs keff'
#ax4.set_title(fig_subtitle)
#ax4.set_xlabel(u'$K_{eff}$'+ u' $(\u00B5M^{-1} \cdot s^{-1})$')
#ax4.set_ylabel('Doubling time (Hrs)')


getScatter (df, ax1, "kcat (1/s)", "Selcoeff" )
getScatter (df, ax2, "kcat/km (uM-1*s-1)", "Selcoeff" )
#getScatter (df, ax3, "kcat (1/s)", "doubling_time" )
#getScatter (df, ax4, "kcat/km (uM-1*s-1)", "doubling_time" )

fig.tight_layout()
fig.subplots_adjust(top = 0.88)
fig.show()


# In[5]:


fig.savefig(fout)

