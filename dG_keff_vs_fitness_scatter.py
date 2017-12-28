
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
df


# In[3]:


a = [["dG_NI", "dG_total"],["kcat (1/s)","kcat/km (10-3 uM-1*s-1)"], ["Selcoeff"]]
xyzcombo = list(itertools.product(*a))


# In[4]:


def ceil_power_of_10(n):
    exp = math.log(abs(n), 10)
    exp = math.ceil(exp)
    return 10**exp

def autolabel(p, x, y, ax, intercept, slope):
    for prot, x_val, y_val in zip(p, x, y):
        ypred = intercept+slope*x_val
        y_shift = ceil_power_of_10(min(y)) * 0.02
        if (y_val >= ypred):
            ax.text(x_val, y_val+y_shift, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val-y_shift, prot, ha='center', va='top')
    return


# In[5]:


def getScatter (df, plotnum, cols, labels, ax_limits): 
    x = df[cols[0]]
    y = df[cols[1]]
    z = df[cols[2]]
    ax = fig.add_subplot(2, 2, plotnum, projection='3d')
    for i,j,k, c in zip(x,y,z, df["color_p"]):
        ax.plot([i,i],[j,j],[k,0], color = c)
    ax.scatter(x, y, c= df["color_p"] , marker = ".", zdir='z', zs = ax_limits[2][0], alpha = 0.3)
    ax.scatter(y, z, c= df["color_p"] , marker = ".", zdir='x', zs = ax_limits[0][0], alpha = 0.3)
    ax.scatter(x, z, c= df["color_p"] , marker = ".", zdir='y', zs = ax_limits[1][1], alpha = 0.3)

    ax.scatter(x, y, z, zdir='z', c= df["color_p"] , marker="o", s=30, linewidth='0.5', edgecolor='black', depthshade = False)
    fig_subtitle = labels[0] + " " + labels[1] + " " +labels[2]
    #x.set_title(fig_subtitle, fontsize=8)
    ax.set_xlabel(labels[0], labelpad = 6)
    ax.set_ylabel(labels[1], labelpad = 6)
    ax.set_zlabel(labels[2], labelpad = 14)
    #ax.grid(False)
    ax.xaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "lightgrey", "lightstyle": ":"})
    ax.yaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "lightgrey", "lightstyle": ":"})
    ax.zaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "lightgrey", "lightstyle": ":"})
    ax.xaxis.pane.set_edgecolor('white')
    ax.yaxis.pane.set_edgecolor('white')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = "whitesmoke"
    #						axes min, max, step
    ax.set_xticks(np.arange(ax_limits[0][0],ax_limits[0][1], ax_limits[0][2]))
    ax.set_yticks(np.arange(ax_limits[1][0],ax_limits[1][1], ax_limits[1][2]))
    ax.set_zticks(np.arange(ax_limits[2][0],ax_limits[2][1], ax_limits[2][2]))
    ax.set_xlim(ax_limits[0][0], ax_limits[0][1])
    ax.set_ylim(ax_limits[1][0], ax_limits[1][1])
    ax.set_zlim(ax_limits[2][0], ax_limits[2][1])
    
    ax.tick_params(axis='z', pad=10)
    ax.tick_params(axis='y', pad=0)
    ax.tick_params(axis='x', pad=0)
    
    #ax.zaxis.labelpad = 20
    #ax.yaxis.labelpad = 15
    
    xx, yy = np.meshgrid([ ax_limits[0][0], ax_limits[0][1]], [ax_limits[1][0], ax_limits[1][1]])
    zz = yy*0
    ax.plot_surface(xx, yy, zz, color = "lightgrey", alpha = 0.1)
    #ax.view_init(elev=20., azim=35)
    ax.yaxis.set_rotate_label(90) 
    #print ax.get_zmajorticklabels()
    xtick = ax.get_xticks()
    ytick = ax.get_yticks()
    ztick = ax.get_zticks()

    ax.plot_wireframe(xx, yy, zz, color = "grey", alpha = 0.2)


    return 


# In[6]:


fig= plt.figure(figsize=(14, 14))
#fig2= plt.figure(figsize=(12, 12))
#fig3 = plt.figure(figsize=(12, 12))

for i, xyz in enumerate(xyzcombo):
    x_col = str(xyz[0])
    y_col = str(xyz[1])
    z_col = str(xyz[2])
    x_label = labeldict[x_col]
    y_label = labeldict[y_col]
    z_label = labeldict[z_col]
    x_ax = axesdict[x_col]
    y_ax = axesdict[y_col]
    z_ax = axesdict[z_col]
    cols = [x_col, y_col, z_col]
    labels = [x_label, y_label, z_label]
    ax_limits = [x_ax, y_ax, z_ax]
    getScatter (df, i+1, cols, labels, ax_limits)

fig_title = u'Stability vs fitness vs enzyme activity' 
fig.suptitle(fig_title, fontsize=12)
plt.subplots_adjust(hspace=0.1)

plt.tight_layout()
fig.show()


# In[7]:


fout = "output/keff_vs_stability_vs_selecoeff"
fig.savefig(fout+".pdf")

