
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


fout = 'output/multiple_regression_dGni_keff.pdf'

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

df["kcat/km (10-3 uM-1*s-1)"] = df["kcat/km (uM-1*s-1)"] * 1e3

df.rename(columns={'dG_NI': 'dg', 'kcat/km (10-3 uM-1*s-1)': 'keff'}, inplace=True)


# In[3]:


df["dgsq"] = df["dg"]**2

X = df[["keff", "dg", "dgsq"]]
y = df["Selcoeff"]
## fit a OLS model 
X = sm.add_constant(X)
est = sm.OLS(y, X).fit()
est.summary()


# In[4]:


est = smf.ols(formula='Selcoeff ~ keff + dg + dgsq', data=df).fit()
print est.params[2], est.params[3]


# In[5]:


from operator import add
#map(add, list1, list2)

## Create the 3d plot --
# dg/keff grid for 3d plot
xx1, xx2 = np.meshgrid(np.linspace(X.keff.min(), X.keff.max(), 100), 
                       np.linspace(X.dg.min(), X.dg.max(), 100))

# calculate quadratic back into single x 
xquad = map(add, est.params[2] * xx2, est.params[3] * xx2) 
#xquad1 = list(chain(*xquad))
#print est.params[1] * xx1
#xquadlist = np.array(xquad.tolist())
# plot the hyperplane by evaluating the parameters on the grid
Z = est.params[0] + est.params[1] * xx1 + xquad


# In[6]:


#The plot above shows data points above the hyperplane in white and points below the hyperplane in black. 
#The color of the plane is determined by the corresonding predicted fitness values (blue = low, red = high). 

fig = plt.figure(figsize=(12, 8))
ax = Axes3D(fig, azim=-115, elev=15)
surf = ax.plot_surface(xx1, xx2, Z, cmap=plt.cm.RdBu_r, alpha=0.6, linewidth=0)
resid = y - est.predict(X)
ax.scatter(X[resid >= 0].keff, X[resid >= 0].dg, y[resid >= 0], color='black', alpha=1.0, facecolor='white')
ax.scatter(X[resid < 0].keff, X[resid < 0].dg, y[resid < 0], color='black', alpha=1.0)
ax.set_ylabel(u'$\u0394G_{NI}$', labelpad = 12)
ax.set_xlabel(u'K$_{eff}$'+ u' $(10^{-3} \cdot \u00B5M^{-1} \cdot s^{-1})$', labelpad = 16)
ax.set_zlabel('s')


# In[7]:


fig.savefig(fout)

