#!/usr/bin/env ipython

# -*- coding: utf-8 -*-

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
from matplotlib.offsetbox import AnchoredText
import itertools
import collections 
import matplotlib as mpl
import matplotlib.artist as artists
from scipy import stats

#get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
from figure_dict import *





#input file - thermal melt first derivative	
fin = "output/Thermalmelt-Tm.csv"
Tm = pd.read_csv(fin)

fin = "output/Titration-ddG.csv"
dG = pd.read_csv(fin)
dG["color_p"] = dG["Protein"].map(colordict)
dG["label_p"] = dG["Protein"].str.replace('_','\n')


#Output
dirOut = "output/"
analysis = "dG_vs_Tm"
fout = dirOut+analysis+'-scatter.pdf'

df = Tm.set_index('Protein').join(dG.set_index('Protein'), how = "inner")
df = df.sort_values(by=['Tm'])



def autolabel(p, x, y, ax, intercept, slope):
    for prot, x_val, y_val in zip(p, x, y):
        ypred = intercept+slope*x_val
        if (y_val >= ypred):
            ax.text(x_val, y_val+0.1, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val-0.1, prot, ha='center', va='top')			
    return

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
    ax.tick_params(axis='both') 

    #ax.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black', label='_nolegend_')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    ax.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black')
    ax.plot(x, intercept+slope*x, '--k')
    ax.text(min(x), max(y), u'$R= %.2f$' % r_value, ha='left', va='center', color = "black")
    
    autolabel(lab_p, x, y, ax, intercept, slope)
    ax.set_xlim([min(x)-1,max(x)+1])
    ax.set_ylim([min(y)-1,max(y)+1])
    return





x_col = "Tm"
y_list = ["dG_NI", "dG_IU", "dG_total"]


fig = plt.figure(figsize=(16, 16)) 
gs = gridspec.GridSpec(nrows=2, ncols=2)
fig_title = u'dG vs Tm' 
fig.suptitle(fig_title)
    
for i, y_col in enumerate (y_list):
    x_label = labeldict[x_col]
    y_label = labeldict[y_col]
    ax = plt.subplot(gs[i])
    fig_subtitle = y_label + " vs " + x_label
    ax.set_title(fig_subtitle)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    scatter(df, ax, x_col, y_col)


fig.tight_layout()
plt.subplots_adjust(hspace=0.25)
fig.subplots_adjust(top = .92)

fig.savefig(fout)

