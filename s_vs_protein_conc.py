#!/usr/bin/env ipython

# coding: utf-8



import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from figure_dict import *


plt.style.use('figures')





#input
dirIn = "input/GrowthAssay/"
fname = "SsMutant"
extracaption = "-unlabeled"
temp = 30 
fin =  dirIn+fname+".xlsx"
xls = pd.ExcelFile(fin)
df = xls.parse("western")

#edit df
df["color_p"] = df["Protein"].map(colordict)
df["label_p"] = df["Protein"].str.replace('_','\n')
#drop vector and any other rows with missing values 
df = df.dropna()
df = df.set_index('Protein')

#output
dirOut = "output/"
analysis = 's_vs_protein_conc'
fout = dirOut+analysis+".pdf"







def autolabel(p, x, y, ax):
    for prot, x_val, y_val in zip(p, x, y):
        invertlabel = ["SsWT", "E155A"]
        moveleftlabel = ["I45K"]
        moverightlabel = ["V78A"]
        moverightuplabel = ["D165A"]
        if (prot in invertlabel): 
            ax.text(x_val, y_val-0.02, prot, ha='center', va='top')
        elif (prot in moveleftlabel): 
            ax.text(x_val-0.05, y_val+0.01, prot, ha='center', va='bottom')
        elif (prot in moverightlabel): 
            ax.text(x_val+0.1, y_val-0.03, prot, ha='center', va='bottom')
        elif (prot in moverightuplabel): 
            ax.text(x_val+0.1, y_val+0.01, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val+0.01, prot, ha='center', va='bottom')
    return

def scatter(ax, x, y, col_p, p, x_col, y_col):
    
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
    
    autolabel(p, x, y, ax)
    padding = 0.1
    ax.set_xlim([min(x)-padding,max(x)+padding])
    ax.set_ylim([min(y)-padding,max(y)+padding])
    #plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    ax.tick_params(axis='both') 
    return





x_col = "ProteinExpression"
y_col = "Selcoeff_all"

x_label = labeldict[x_col]
y_label = labeldict[y_col]

fig = plt.figure(figsize=(8,8)) 
ax = plt.subplot()
fig_title = y_label + " vs " + x_label
fig.suptitle(fig_title)
plt.subplots_adjust(hspace=0.1)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)

x = df[x_col]
y = df[y_col]
col_p = df["color_p"]
p = df["label_p"]

scatter(ax, x, y, col_p, p, x_col, y_col)

fig.tight_layout()
fig.subplots_adjust(top = 0.88)



fig.savefig(fout)

