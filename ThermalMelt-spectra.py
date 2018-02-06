#!/usr/bin/env ipython
# coding: utf-8



import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
import itertools
import collections 
import matplotlib as mpl
import matplotlib.artist as artists
#get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
from figure_dict import *


#input file - thermal melt first derivative	
fin = "input/ThermalMelt/20170725-Derivative.xlsx"
xls = pd.ExcelFile(fin)
df = xls.parse("norm_deriv")
extracaption = "TexasRed-FirstDerivative"

#Output
dirOut = "output/"
experiment = "Thermalmelt"
fout = dirOut+experiment+'-spectra.pdf'
TMout = dirOut+experiment+'-Tm.csv'


df.head(3)





#save temp to use as x-values 
x = df["Temp"]
#remove temp colun from df 
df =   df.iloc[:,1:]

dict_temp1 = df.to_dict(orient = "list")
protlist = list(df)





#read in list of files, use to get proteins  
marker_p =   ['o'] * len(protlist)
markersize_p =  [20] * len(protlist)
markerdict = dict(zip(protlist, marker_p))
sizedict = dict(zip(protlist, markersize_p))

order={}
for p in protlist: 	
	order[orderdict[p]] = p

sorted_protlist = []

for key, value in order.iteritems():	
	sorted_protlist.append(value) 




###--- plot scatter and line graph of urea titration 
def autolabel(colorprot, tm, maxvalue):
    height = maxvalue+0.001
    if colorprot == "yellow":
        colorprot = "black"
    ax1.text(x = tm, y = height+0.002, s = str(tm), ha='center', va='center', color = colorprot, fontsize=16 )





fig = plt.figure(figsize=(8,8)) 
ax1 = plt.subplot()

fig_title = u'Thermal Melt - SYPRO RED, Ex/EM Texas Red       ' 
fig.suptitle(fig_title)
plt.subplots_adjust(hspace=0.1)
fig_subtitle = extracaption

ax1.set_xlabel(labeldict["Temp"])
ax1.set_ylabel(labeldict["RFU-Temp"])

tm_dict = {}
for p in sorted_protlist: 
    ax1.plot(x, dict_temp1[p], ls = "-", lw = 1,  c= colordict[p], label='_nolegend_')
    ax1.scatter(x, dict_temp1[p], c= colordict[p], marker=markerdict[p], s=sizedict[p], 
                label=p.replace('_', '/'), linewidth='0.5', edgecolor='black')
    maxvalue = max(dict_temp1[p])
    #tm = max(enumerate(dict_temp1[p]), key=lambda x: x[1])[0]
    [tm for y, tm in sorted(zip(dict_temp1[p], x))] 
    tm_dict[p] = tm
    #max([y for x,dict_temp1[p] in ziplist])
    autolabel(colordict[p], tm, maxvalue)
ax1.tick_params(axis='both') 
##, bbox_to_anchor=(1, 1)
ax1.legend(scatterpoints=1, loc="best", handletextpad=0.05, fontsize=14, frameon=False)
ax1.set_xlim([x.min(),x.max()])	
fig.tight_layout()
fig.subplots_adjust(top = 0.88)
#fig.show()





fig.savefig(fout)




#save tm to file 
with open( TMout, 'wb') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["Protein", "Tm"])
    for key, value in tm_dict.items():
       writer.writerow([key, value])

