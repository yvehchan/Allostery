#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib.ticker as ticker
import itertools
import collections 
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties import unumpy
from figure_dict import *
import os

plt.style.use('figures')


# temp = 30C
# experiment = titration 


#Input
dirIn = "input/CD-titration/"

#description 
fname = "SsMutant"
temp =30

#Output
dirOut = "output/"
experiment = "Titration"
fout = dirOut+experiment+'-spectra.pdf'
ddGout = dirOut+experiment+'-ddG.csv'
addout = dirOut+experiment+'-additivity.csv'
combinedout = dirOut+experiment+'-ddG-additivity.csv'
dGplotout = dirOut+experiment+'-bargraph-dG.pdf'
ddGplotout = dirOut+experiment+'-bargraph-ddG.pdf'
addplotout = dirOut+experiment+'-bargraph-additivity.pdf'


urea = {}
MREdata = {}
urea_model = {}
model = {}
urea_residual = {}
residual = {}
gni = {}
gni_err = {}
giu = {}
giu_err = {}
gtotal = {} 

#read in list of files, use to get proteins  
protlist = []
minvalue = 10000
maxvalue = -10000


# Inputs are CD spectra and their fits determined using Savuka


filelist = [x for x in os.listdir(dirIn) if x.endswith(".xlsx")]
#print(filelist)




for line in filelist:
    line = line.strip()
    _, p, _, _ = line.split('-')
    if (p != "SsWT"):
        p = p[2:]
    #print (p)
    protlist.append(p)
    xls = pd.ExcelFile(dirIn+line)
    df = xls.parse('fit') #titration values and fit
    urea[p] = df.iloc[:,0].dropna().values
    #temp = np.divide(df.iloc[:,2].dropna().values,1000)
    MREdata[p] = np.divide(df.iloc[:,2].dropna().values,1000)
    urea_model[p] = df.iloc[:,4].dropna().values
    model[p] = np.divide(df.iloc[:,6].dropna().values,1000)
    urea_residual[p] = df.iloc[:,8].dropna().values
    residual[p] = np.divide(df.iloc[:,9].dropna().values,1000)
    tempmin = min(np.divide(df.iloc[:,9].dropna().values,1000))
    tempmax = max(np.divide(df.iloc[:,9].dropna().values,1000))
    if  tempmin < minvalue:
        minvalue = tempmin
    if  tempmax > maxvalue:
        maxvalue = tempmax
    df2 = xls.parse('dg', header = None) #dg only
    gni[p] = df2.iloc[0,1] #row, column
    gni_err[p] = df2.iloc[0,2] 
    giu[p] = df2.iloc[1,1] 
    giu_err[p] = df2.iloc[1,2] 
    gtotal[p] = df2.iloc[0,1] + df2.iloc[1,1] 




marker_p =   ["o"] * len(protlist)
markersize_p =  [50] * len(protlist)
markerdict = dict(zip(protlist, marker_p))
sizedict = dict(zip(protlist, markersize_p))


order={}
for p in protlist:
	order[orderdict[p]] = p

sorted_protlist = []
sorted_colorlist = []
for key, value in order.items():
	sorted_protlist.append(value) 




urea["SsWT"].size


fig = plt.figure(figsize=(12, 12)) 
gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[4,1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
fig_title = u'Urea titration - MRE vs Urea [M]  ' 
fig.suptitle(fig_title)
plt.subplots_adjust(hspace=5)
fig_subtitle = str(temp) + '$^\circ$C' + ", 10mM KPi pH 7.2"
ax1.set_title(fig_subtitle)
ax1.set_xlabel(u'[Urea] (M)',)
#ax1.set_ylabel(u'\$MRE ^222 ^nm$(deg*cm$\_2*dmol\_-1$' \)', fontsize=16)
ax1.set_ylabel(r'$MRE_{222}(10^3 deg \cdot cm^2 \cdot dmol^{-1})$')
ax2.set_xlabel(u'[Urea] (M)')
ax2.set_ylabel('Residual')

for p in sorted_protlist: 	
    if p == "SsWT":
        ax1.plot(urea_model[p], model[p], ls = "--", lw = 0.5,  c= "black", label='_nolegend_')
        ax1.scatter(urea[p], MREdata[p], c= "white", marker=markerdict[p], s=sizedict[p], label=p, linewidth='1', edgecolor='black')
    else:
        ax1.plot(urea_model[p], model[p], ls = "--", lw = 0.5,  c= colordict[p], label='_nolegend_')
        ax1.scatter(urea[p], MREdata[p], c= colordict[p], marker=markerdict[p], s=sizedict[p], label=p.replace('_', '/'), linewidth='0.5', edgecolor='black')
    ax2.scatter(urea_residual[p], residual[p], c= colordict[p], marker=markerdict[p], s=sizedict[p], linewidth='0.5', edgecolor='black')
#ax1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
ax1.legend(loc="lower right", scatterpoints = 1, frameon=False,handletextpad=0.05)
ax1.set_xlim([0,10])	
ax2.set_xlim([0,10])	
ax2.set_ylim([minvalue-0.001,maxvalue+0.001])	
fig.tight_layout()
fig.subplots_adjust(top = 0.88, right=0.8)
#fig.show()





#print dirOut+fout
fig.savefig(fout)




def propagateErrorDF(unumpyarray, df, colname):
    unumpyarray = ['{:.6f}'.format(i) for i in unumpyarray]
    tempSeries = pd.Series(unumpyarray, name = colname).apply(str)
    tempDF = pd.DataFrame([ x.split('+/-') for x in tempSeries.tolist() ], columns = [colname, colname + '_err'], index = df.index)
    df = df.join(tempDF) 
    return (df)





gni_list_values = []
gni_err_list_values = []
giu_list_values = []
giu_err_list_values = []
gtotal_list_values = []
sorted_colorlist = []




for p in sorted_protlist:
    gni_list_values.append(gni[p])
    gni_err_list_values.append(gni_err[p])
    giu_list_values.append(giu[p])
    giu_err_list_values.append(giu_err[p])
    gtotal_list_values.append(gtotal[p])
    sorted_colorlist.append(colordict[p])


# #save as ifloat to propagate uncertainties, variables ending with "_u" = uncertainties format (value, err)
# #ddG = dGmut-dGwt, more neg dG is more stable, ddG >0 = less stable than reference
# #dGs are neg, but shown without neg sign, so need to multiple values by -1 
# #substract mut from wt for equivalent, corrected sign 




ddG = pd.DataFrame(np.column_stack([sorted_protlist, gni_list_values,gni_err_list_values,giu_list_values,giu_err_list_values,gtotal_list_values]), columns=['Protein', 'dG_NI', 'dG_NI_err', 'dG_IU', 'dG_IU_err', 'dG_total'])
ddG = ddG.set_index("Protein")

ddG_NI_u = ufloat( ddG["dG_NI"]["SsWT"], ddG["dG_NI_err"]["SsWT"]) - unumpy.uarray(ddG["dG_NI"].values, ddG["dG_NI_err"].values)
ddG_IU_u = ufloat( ddG["dG_IU"]["SsWT"], ddG["dG_IU_err"]["SsWT"]) - unumpy.uarray(ddG["dG_IU"].values, ddG["dG_IU_err"].values)
#ddG_NI_u = unumpy.uarray(ddG["dG_NI"].values, ddG["dG_NI_err"].values) - ufloat( ddG["dG_NI"]["SsWT"], ddG["dG_NI_err"]["SsWT"])
#ddG_IU_u = unumpy.uarray(ddG["dG_IU"].values, ddG["dG_IU_err"].values) - ufloat( ddG["dG_IU"]["SsWT"], ddG["dG_IU_err"]["SsWT"])

 
ddG = propagateErrorDF(ddG_NI_u, ddG, "ddG_NI")
ddG = propagateErrorDF(ddG_IU_u, ddG, "ddG_IU")

ddG = ddG.astype('float')	
ddG["ddG_total"] = ddG["dG_total"]["SsWT"] - ddG["dG_total"]
ddG.to_csv(ddGout, index=True)





def getInteractionEnergy(p):
    mut1, mut2 = p.split("_")
    ## ufloat( value, err )
    double_ni_u = ufloat(ddG["ddG_NI"][p], ddG["ddG_NI_err"][p])
    mut1_ni_u   = ufloat(ddG["ddG_NI"][mut1], ddG["ddG_NI_err"][mut1])
    mut2_ni_u   = ufloat(ddG["ddG_NI"][mut2], ddG["ddG_NI_err"][mut2])
    #print mut1, mut2, mut1_ni_u, mut2_ni_u, double_ni_u
    deltaNI_u = -double_ni_u - (-mut1_ni_u + -mut2_ni_u)
    deltaNI_u = '{:.6f}'.format(deltaNI_u)
    deltaNI, deltaNI_err = deltaNI_u.split("+/-")
    double_iu_u = ufloat(ddG["ddG_IU"][p], ddG["ddG_IU_err"][p])
    mut1_iu_u   = ufloat(ddG["ddG_IU"][mut1], ddG["ddG_IU_err"][mut1])
    mut2_iu_u   = ufloat(ddG["ddG_IU"][mut2], ddG["ddG_IU_err"][mut2])
    deltaIU_u = -double_iu_u - (-mut1_iu_u + -mut2_iu_u)
    deltaIU_u = '{:.6f}'.format(deltaIU_u)
    deltaIU, deltaIU_err = deltaIU_u.split("+/-")	
    deltaTotal = ddG["ddG_total"][p] - (ddG["ddG_total"][mut1] + ddG["ddG_total"][mut2])
    NI = [ddG["ddG_NI"][mut1], ddG["ddG_NI_err"][mut1], ddG["ddG_NI"][mut2], ddG["ddG_NI_err"][mut2], ddG["ddG_NI"][p], ddG["ddG_NI_err"][p], deltaNI, deltaNI_err]
    IU = [ddG["ddG_IU"][mut1], ddG["ddG_IU_err"][mut1], ddG["ddG_IU"][mut2], ddG["ddG_IU_err"][mut2], ddG["ddG_IU"][p], ddG["ddG_IU_err"][p], deltaIU, deltaIU_err]
    Total = [ddG["ddG_total"][mut1], ddG["ddG_total"][mut2], ddG["ddG_total"][p], deltaTotal] 
    results = [p] + NI + IU + Total
    return (results)




additivity_list = ["I45A_S70A", "I45A_M73A"]
additivity = pd.DataFrame(columns=["Protein", 
									"ddG_NI_mut1", "ddG_NI_err_mut1", "ddG_NI_mut2", "ddG_NI_err_mut2", "ddG_NI_double", "ddG_NI_err_double","interaction_NI", "interaction_NI_err",
									"ddG_IU_mut1", "ddG_IU_err_mut1", "ddG_IU_mut2", "ddG_IU_err_mut2", "ddG_IU_double", "ddG_IU_err_double","interaction_IU", "interaction_IU_err",
									"ddG_total_mut1", "ddG_total_mut2","ddG_total_double", "interaction_total" ])
for p in additivity_list:
    results = getInteractionEnergy(p)
    additivity.loc[additivity.shape[0]] = results
additivity = additivity.set_index("Protein")
#print additivity
additivity.to_csv(addout, index=True)



combinedTable = pd.DataFrame() 
combinedTable["dG_NI_plus_err"] = ddG["dG_NI"].round(2).map(str) + u'\u00B1' + ddG["dG_NI_err"].round(2).map(str)
combinedTable["dG_IU_plus_err"] = ddG["dG_IU"].round(2).map(str) + u'\u00B1' + ddG["dG_IU_err"].round(2).map(str)
combinedTable["ddG_NI_plus_err"] = ddG["ddG_NI"].round(2).map(str) + u'\u00B1' + ddG["ddG_NI_err"].round(2).map(str)
combinedTable["ddG_IU_plus_err"] = ddG["ddG_IU"].round(2).map(str) + u'\u00B1' + ddG["ddG_IU_err"].round(2).map(str)
combinedTable["interaction_NI_plus_err"] = additivity["interaction_NI"].astype(float).round(2).map(str) + u'\u00B1' + additivity["interaction_NI_err"].astype(float).round(2).map(str)
combinedTable["interaction_IU_plus_err"] = additivity["interaction_IU"].astype(float).round(2).map(str) + u'\u00B1' + additivity["interaction_IU_err"].astype(float).round(2).map(str)
combinedTable = combinedTable.reindex(sorted(combinedTable.index, key=lambda x: orderdict[x]))
combinedTable['ddG_NI_plus_err']['SsWT'] = "-"
combinedTable['ddG_IU_plus_err']['SsWT'] = "-"
combinedTable = combinedTable.fillna('-')
#col1 = labeldict['dG_NI']
#col2 = labeldict['dG_IU']
#col3 = labeldict['ddG_NI']
#col4 = labeldict['ddG_IU']
#col5 = labeldict['interaction_NI']
#col6 = labeldict['interaction_IU']
combinedTable.columns = ["dG_NI", "dG_IU", "ddG_NI", "ddG_IU",  "interaction_NI", "interaction_IU"]

combinedTable.to_csv(combinedout, encoding='utf-8-sig', index=True)
#combinedTable





def autolabel(rects, list_values):
	for rect, h in zip(rects, list_values):
		height = h
		if (height >= 0):
			ax.text(rect.get_x() + rect.get_width()/2, height+.3, '%.2f' % height, 
                    ha='center', va='center', fontsize=16)
		if (height < 0):
			ax.text(rect.get_x() + rect.get_width()/2, height-.3, '%.2f' % height, 
                    ha='center', va='center', fontsize=16)

####---bargraph for dG(ni) and dG(iu)				
#for bargraph change WT color to white 
combined_list_values = gni_list_values + giu_list_values
sorted_colorlist_Gni = ["white" if x=="black" else x for x in sorted_colorlist]
sorted_colorlist_Giu = ["white" if x=="black" else x for x in sorted_colorlist]
N = len(sorted_protlist)
labels = [w.replace('_', '\n') for w in sorted_protlist]
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
fig = plt.figure(figsize=(12, 8)) 
gs = gridspec.GridSpec(nrows=1, ncols=1)
ax = plt.subplot(gs[0])
rects1 = ax.bar(ind, gni_list_values, width, color=sorted_colorlist_Gni, yerr=gni_err_list_values,edgecolor='black')
plt.rcParams['hatch.linewidth'] = 3
rects2 = ax.bar(ind + width, giu_list_values, width, color=sorted_colorlist_Giu, yerr=giu_err_list_values,edgecolor='black', hatch='//')
# add some text for labels, title and axes ticks
ax.set_ylabel(u'\u0394G ' + r'$(kcal \cdot mol^{-1})$')

ax.set_title(u'\u0394G ' + ' folding free energy')
ax.set_xticks(ind + width / 2)
ax.tick_params(axis='both') 
ax.set_xticklabels(labels)
ax.set_ylim([0, max(combined_list_values)*1.2])
ax.legend((rects1[0], rects2[0]), (u'$\u0394G_{NI}$', u'$\u0394G_{IU}$'),  prop={'size': 20})
autolabel(rects1, gni_list_values)
autolabel(rects2, giu_list_values)
ax.axhline(y=0, color = "black")
#fig.show()
fig.savefig(dGplotout)


####---bargraph for ddG(ni) and ddG(iu)	

ddG = ddG.filter(regex='ddG')
ddG = ddG.drop('SsWT', axis=0)
ddG["Prot"] = ddG.index
ddG['color'] = ddG["Prot"].map(colordict)
ddG['ddG_NI_plus'] = ddG["ddG_NI"]+ddG["ddG_NI_err"]
ddG['ddG_NI_minus'] = ddG["ddG_NI"]-ddG["ddG_NI_err"]
ddG['ddG_IU_plus'] = ddG["ddG_IU"]+ddG["ddG_IU_err"]
ddG['ddG_IU_minus'] = ddG["ddG_IU"]-ddG["ddG_IU_err"]
sorted_protlist = ddG["Prot"].tolist()
labels = [w.replace('_', '\n') for w in sorted_protlist]
sorted_colorlist = ddG['color'].tolist()
sorted_colorlist_Giu = ["white" if x=="black" else x for x in sorted_colorlist]
ddgni_list_values =  ddG['ddG_NI'].tolist()
ddgni_err_list_values =  ddG['ddG_NI_err'].tolist()
ddgiu_list_values =  ddG['ddG_IU'].tolist()
ddgiu_err_list_values =  ddG['ddG_IU_err'].tolist()
## sue to get max/min values of y-axis
combined_list_values =  ddG['ddG_NI_plus'].tolist() + ddG['ddG_NI_minus'].tolist() + ddG['ddG_IU_plus'].tolist() + ddG['ddG_IU_minus'].tolist() 
ddgtotal_list_values =  ddG['ddG_total'].tolist()




N = len(sorted_protlist)
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
fig = plt.figure(figsize=(12, 8)) 
gs = gridspec.GridSpec(nrows=1, ncols=1)
ax = plt.subplot(gs[0])
##dummy plot for legend marker
rects3 = ax.bar(ind, ddgni_list_values, width, color="white",edgecolor='black')
rects4 = ax.bar(ind + width, ddgiu_list_values, width, color="white", edgecolor='black', hatch='//')
rects1 = ax.bar(ind, ddgni_list_values, width, color=sorted_colorlist,  yerr=ddgni_err_list_values, edgecolor='black')
plt.rcParams['hatch.linewidth'] = 3
rects2 = ax.bar(ind + width, ddgiu_list_values, width, color=sorted_colorlist_Giu, yerr=ddgiu_err_list_values, edgecolor='black', hatch='//')
# add some text for labels, title and axes ticks
ax.set_ylabel(u'\u0394\u0394G ' + r'$(kcal \cdot mol^{-1})$', fontsize=20)
ax.set_title(u'\u0394\u0394G ' + r'$(kcal \cdot mol^{-1})$' + ' difference in folding free energy between wild type and mutant')
ax.set_xticks(ind + width / 2)
ax.tick_params(axis='both') 
ax.set_xticklabels(labels)
ax.set_ylim([min(combined_list_values)*1.5, max(combined_list_values)*1.2])
autolabel(rects1, ddgni_list_values)
autolabel(rects2, ddgiu_list_values)
ax.legend((rects3[0], rects4[0]), (u'$\u0394\u0394G_{NI}$', u'$\u0394\u0394G_{IU}$'),  prop={'size': 20}, loc=3)
ax.axhline(y=0, color = "black")
#fig.show()




fig.savefig(ddGplotout)





#fin = "SsMutantTitrations-25C-analysis-additivity.csv"
#add = pd.read_csv(fin)
add = additivity

add['Protein'] = add.index
#add.convert_objects(convert_numeric=True)
add = add.apply(pd.to_numeric, errors="ignore")


ddgni_list_values =  add['interaction_NI'].astype(float).tolist()
ddgni_err_list_values =  add['interaction_NI_err'].astype(float).tolist()
ddgiu_list_values =  add['interaction_IU'].astype(float).tolist()
ddgiu_err_list_values =  add['interaction_IU_err'].astype(float).tolist()




minvalue = min(min(ddgni_list_values), min(ddgiu_list_values))
maxvalue = max(max(ddgni_list_values), max(ddgiu_list_values))
protlist = add['Protein'].tolist()
labels = [w.replace('_', '\n') for w in protlist]

def autolabel(rects, list_values):
	for rect, h in zip(rects, list_values):
		height = h
		if (height >= 0):
			ax.text(rect.get_x() + rect.get_width()/2, height+.3, '%.2f' % height, 
                    ha='center', va='center', fontsize=16)
		if (height < 0):
			ax.text(rect.get_x() + rect.get_width()/2, height-.3, '%.2f' % height, 
                    ha='center', va='center', fontsize=16)




colorlist = ["green", "blue"]
#mpl.rcParams['hatch.linewidth'] = 3

N = 2
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
fig = plt.figure(figsize=(12, 8)) 
#plt.rcParams["patch.force_edgecolor"] = True
gs = gridspec.GridSpec(nrows=1, ncols=1)
ax = plt.subplot(gs[0])

##dummy plot for legend marker
rects3 = ax.bar(ind, ddgni_list_values, width, color="white",edgecolor='black')
rects4 = ax.bar(ind + width, ddgiu_list_values, width, color="white", edgecolor='black', hatch='//')


rects1 = ax.bar(ind, ddgni_list_values, width, color=colorlist, yerr=ddgni_err_list_values, edgecolor='black')
rects2 = ax.bar(ind+width, ddgiu_list_values, width, color=colorlist, yerr=ddgiu_err_list_values)
rects5 = ax.bar(ind+width, ddgiu_list_values, width, color="None", hatch='//')

# add some text for labels, title and axes ticks
ax.set_ylabel(u'\u03B4 ' + r'$(kcal \cdot mol^{-1})$')
ax.tick_params(axis='both') 
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(labels)
#ax.set_ylim([minvalue*1.5, maxvalue*1.2])
autolabel(rects1, ddgni_list_values)
autolabel(rects2, ddgiu_list_values)

ax.legend((rects3[0], rects4[0]), (u'$\u03B4_{NI}$', u'$\u03B4_{IU}$'),  prop={'size': 20}, loc = 4)
ax.axhline(y=0, color = "black")



fig.savefig(addplotout)



