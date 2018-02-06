#!/usr/bin/env ipython

# -*- coding: utf-8 -*-




import csv
import pandas as pd
import numpy as np
import itertools
import collections 
from scipy import stats
import math
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
from matplotlib.offsetbox import AnchoredText
import matplotlib as mpl
import matplotlib.artist as artists
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from figure_dict import *

plt.style.use('figures')

#input
fin = "output/Titration-ddG.csv"
dG = pd.read_csv(fin)

fin = "input/GrowthAssay/SsMutant.xlsx"
xls = pd.ExcelFile(fin)
dbl = xls.parse("selcoeff")

fin = "output/Thermalmelt-Tm.csv"
tm = kcat = pd.read_csv(fin)
tm["color_p"] = tm["Protein"].map(colordict)
tm["label_p"] = tm["Protein"].str.replace('_','\n')

#fin = "output/SsMutant-double-pooled-30C-initialvelocity-mm.csv"
#kcat_temp = "30"
#kcat = pd.read_csv(fin)
#kcat["color_p"] = kcat["Protein"].map(colordict)
#kcat["label_p"] = kcat["Protein"].str.replace('_','\n')

#output
dirOut = "output/"
analysis = 's_vs_stability_dG_Tm'
fout = dirOut+analysis+".pdf"

#joint dataframe 
df = dbl.set_index('Protein').join(dG.set_index('Protein'), how = "inner")
#df = df.join(kcat.set_index('Protein'), how = "inner")
df = df.join(tm.set_index('Protein'), how = "inner")
df

def autolabel(p, x, y, ax, intercept, slope, minx, maxx):
    for prot, x_val, y_val in zip(p, x, y):
		ypred = intercept+slope*x_val
		#ax_ymin, ax_ymax = ax.get_ylim()
		ax_xmin, ax_xmax = ax.get_xlim()
		x_shift = 0.10
		if x_val+x_shift > ax_xmax:
			x_val = x_val - x_shift
		if x_val-x_shift < ax_xmin:
			x_val = x_val + x_shift
		y_shift = 0.005
		if (y_val > ypred):
			ax.text(x_val, y_val+y_shift, prot, ha='center', va='bottom')
		else:
			ax.text(x_val, y_val-y_shift, prot, ha='center', va='top')
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
	ax.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black', label='_nolegend_')

	X = np.vander(x, 2)
	model1 = sm.OLS(y,X)
	results1 = sm.OLS(y,X).fit()
	xx = np.linspace(min(x), max(x), 100)
	XX = np.vander(xx, 2)
	ax.plot(xx, results1.predict(XX), 'k--', label='Poly n=1, $R$=%.2f' % math.sqrt( results1.rsquared), 
			alpha=0.9)

	# 2-nd order polynomial
	# add a constant and quadratic term
	X2 = np.vander(x, 3)
	model2 = sm.OLS(y,X2)
	results2 = sm.OLS(y,X2).fit()
	XX = np.vander(xx, 3)
	ax.plot(xx, results2.predict(XX), 'g--', label='Poly n=2, $R$=%.2f' % math.sqrt( results2.rsquared), 
			alpha=0.9)

	slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
	autolabel(lab_p, x, y, ax, intercept, slope, min(x), max(x))
	if x_col == "dG_NI" or x_col == "dG_total":
		paddingx = 0.4
	elif x_col == "dG_IU":
		paddingx = 0.2
	else:
		paddingx = 1.5
	paddingy = 0.05
	ax.set_xlim([min(x)-paddingx,max(x)+paddingx])
	ax.set_ylim([min(y)-paddingy,max(y)+paddingy])
	ax.tick_params(axis='both') 
	ax.legend(loc = 4)
	return




x_list = ["dG_NI", "dG_IU", "dG_total", "Tm"]
y_col = "Selcoeff"

fig = plt.figure(figsize=(16, 16)) 
gs = gridspec.GridSpec(nrows=2, ncols=2)
fig_title = u'fitness vs stability' 
fig.suptitle(fig_title)
    
for i, x_col in enumerate (x_list):
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

