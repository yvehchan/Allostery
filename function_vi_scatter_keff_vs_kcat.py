
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
import seaborn as sns
from figure_dict import *
import statsmodels.api as sm
import lmfit as lm
import re
from scipy import stats
from scipy.optimize import curve_fit
from collections import defaultdict
import numbers
from uncertainties import ufloat
from uncertainties import umath
from uncertainties.umath import *
from uncertainties import unumpy
from decimal import Decimal

get_ipython().magic(u'matplotlib inline')
plt.style.use('figures')
####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib
##print(plt.style.available)


# In[2]:


def concatDF (dfs):
    df = pd.concat(dfs, join = 'outer', axis =1)
    df.columns = df.columns.droplevel(0)
    df = df.groupby(by=df.columns, axis=1).apply(lambda g: g.mean(axis=1) if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[:,0])  
    return df


# In[3]:


def min_mm(params, x, y):
    V = params["V"].value
    K = params["K"].value
    r = params["r"].value
    model = (V * x/(K + x)) - r
    return model - y

def get_mm(x, V, K):
    Y = (V * x/(K + x))
    return Y

def mm(x, y, V, K, r):
    model = (V * x/(K + x)) - r
    return model 


# In[4]:


def getErrorFromSplit (uncertain_string):
    uncertain_string = '{:.6f}'.format(uncertain_string)
    value, err = uncertain_string.split("+/-")
    return ('%.6E' % Decimal(float(err)))

def getValueFromSplit (uncertain_string):
    uncertain_string = '{:.6f}'.format(uncertain_string)
    value, err = uncertain_string.split("+/-")
    #value = '{:.6ue}'.format(float(value))
    return ('%.6E' % Decimal(float(value)))

def propagateErrorDF(unumpyarray, df, colname):
    unumpyarray = ['{:.20f}'.format(i) for i in unumpyarray]
    tempSeries = pd.Series(unumpyarray, name = colname).apply(str)
    tempDF = pd.DataFrame([ x.split('+/-') for x in tempSeries.tolist() ], columns = [colname, colname + '_err'], index = df.index)
    df = df.join(tempDF)
    return (df)

def autolabel(p, x, y, ax, intercept, slope):
    for prot, x_val, y_val in zip(p, x, y):
        ypred = intercept+slope*x_val
        if (y_val >= ypred):
            ax.text(x_val, y_val+0.00075, prot, ha='center', va='bottom')
        else: 
            ax.text(x_val, y_val-0.00075, prot, ha='center', va='top')
    return


# In[5]:


dirIn = "input/Function/"
dirOut = 'output/'
filelist = ["20170628-SsHelix-3reps", "20170717-SsHelix-plate1-3reps", "20170717-SsHelix-plate2-3reps"]
temp = 30
extracaption = "-"+str(temp)+"C"

vis = {}
stdevs = {}
for i, fname in enumerate(filelist): 
	curr_vi = dirIn+fname+'-initialvelocity-vi.csv'
	vi = pd.read_csv(curr_vi)
	vi.set_index('CdRP[uM]', inplace=True)
	vis[i] = vi 
	curr_std = dirIn+fname+'-initialvelocity-stdev.csv'
	std = pd.read_csv(curr_std)
	std.set_index('CdRP[uM]', inplace=True)
	stdevs[i] = std
velocity = concatDF(vis)	
velstd = concatDF(stdevs)	


# In[6]:


def createScatterWithLinearRegression(df, fout):
    print df
    x = df["kcat (1/s)"].astype('float')
    y = df["kcat/km (uM-1*s-1)"].astype('float')
    col_p = df["color_p"]
    p = df["label_p"]

    fig = plt.figure(figsize=(8, 8)) 
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax1 = plt.subplot(gs[0])
    fig_title = u'$K_{eff}$'+ " vs " + u'$K_{cat}$' 
    #fig.suptitle(fig_title, fontsize=20)
    plt.subplots_adjust(hspace=0.1)
    ax1.set_xlabel(u'$k_{cat}$'+ u' $(s^{-1})$')
    ax1.set_ylabel(u'$k_{eff}$'+ u' $(\u00B5M^{-1} \cdot s^{-1})$')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    ax1.scatter(x, y, c= col_p , marker="o", s=30, linewidth='0.5', edgecolor='black')
    ax1.plot(x, intercept+slope*x, '--k')
    autolabel(p, x, y, ax1, intercept, slope)
    ax1.text(min(x), max(y), u'$R= %.2f$' % r_value, ha='left', va='center')
    ax1.set_xlim([min(x)-1,max(x)+1])	
    ax1.set_ylim([min(y)-0.005,max(y)+0.005])
    wtx = x[0]
    wty = y[0]
    ax1.axvline(x=wtx, color = "grey", linestyle = ":", zorder=0)
    ax1.axhline(y=wty, color = "grey", linestyle = ":", zorder=0)

    fig.subplots_adjust(top = 0.88)
    ax1.tick_params(axis='both') 
    fig.tight_layout()
    fig.show()
    fig.savefig(fout)
    return


# In[7]:


def createScatterErrorBarByGroups(velocity, velstd, group, protlist, proteinconc, extendx = 100):
    velocity.reset_index(level=0, inplace=True)
    velstd.reset_index(level=0, inplace=True)

    fig = plt.figure(figsize=(12, 12)) 
    #size ratio of subplots: table, data, residual from fit
    #gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[2, 4,1])
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5,1])
    #ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    fig_title =  ', '.join(filelist)
    sub_title = group+'-vi-scatter-stdev-residual'
    my_suptitle = plt.suptitle(fig_title)
    #plt.subplots_adjust(hspace=5)
    #ax1.set_title(sub_title, fontsize=18)
    ax1.set_xlabel(u'[CdRP] (\u03bcM)')
    ax1.set_ylabel(u'Initial velocity (\u0394RFU/s)')
    ax2.set_xlabel(u'[CdRP] (\u03bcM)')
    ax2.set_ylabel('Residual')

    x_conc = pd.Series(velocity["CdRP[uM]"])
    ymax = -1

    #get fit of proteins
    resid = dict()
    resid = {"CdRP[uM]": x_conc.tolist()}
    fit = pd.DataFrame(columns = ("Protein", "[Protein] (uM)", "Vmax (uM/s)", "Vmax_err","kcat (1/s)", "kcat_err","Km (uM)", "Km_err","kcat/km (uM-1*s-1)", "kcat/km_err"))

    for p in protlist: 	
        y = pd.to_numeric(velocity[p])
        x = pd.to_numeric(x_conc)
        e = pd.to_numeric(velstd[p])
        xmin = int(x_conc.min())
        xmax = int(x.max()) 
        #concentration to add for extrapolation 
        extended_xmax = xmax + extendx
        #individual outlier removal
        idx = y.index[y.apply(np.isnan)]
        if idx.size:
            for i in idx:
                x = x.drop([i])
                y = y.drop([i])	
                e = e.drop([i])
        #create parameters object
        params = lm.Parameters()
        #initial estimates
        vmax = float(y.max())
        vmin = float(y.min())
        if (vmax > ymax):
            ymax = vmax
        vmaxhalf = (vmax+vmin)/2
        #list of y values close to cmaxhalf
        v_halfmax = y <= vmaxhalf
        #get index of y at half max
        vel_halfmax = y.loc[v_halfmax].idxmax()
        conc_halfmax = float(x.loc[vel_halfmax])
        params.add("V", value = vmax+0.1, max = vmax+1000.0) # vmax estimate
        params.add("K", value = conc_halfmax+0.1, max = conc_halfmax+1000.0)# # km estimate
        params.add("r", value = vmin)
        #Minimize fit of MM, leastsq’: Levenberg-Marquardt (default), 
        result = lm.minimize(min_mm, params, args = (x, y))
        #errorofit = np.sqrt(np.diag(result.covar))
        new_vmax = result.params["V"].value
        new_km = result.params["K"].value
        #standard error = 1 sigma 
        new_vmax_err = result.params["V"].stderr
        new_km_err = result.params["K"].stderr	
        #save as ifloat to propagate uncertainties, variables ending with "_u" = uncertainties format (value, err)
        new_vmax_u = ufloat(new_vmax, new_vmax_err)
        new_km_u = ufloat(new_km, new_km_err)
        kcat_u = new_vmax_u/proteinconc
        keff_u = (new_vmax_u/proteinconc)/new_km_u
        kcat = getValueFromSplit(kcat_u)
        keff = getValueFromSplit(keff_u)
        fit.loc[len(fit)]=(p, proteinconc, round(new_vmax,3), getErrorFromSplit(new_vmax_u), kcat, getErrorFromSplit(kcat_u), round(new_km, 3), getErrorFromSplit(new_km_u),keff, getErrorFromSplit(keff_u))		
        resid.update({p: result.residual})
        extrapolated_x = x.append( pd.Series([extended_xmax]), ignore_index=True)
        extrapolated_y = get_mm(extrapolated_x, new_vmax, new_km)
        #extrapolated_y = y.append( pd.Series([(new_vmax * extrapolated_x/(new_km + extrapolated_x))]), ignore_index=True)
        marker_p =   ["o"] * len(protlist)
        markersize_p =  [20] * len(protlist) 
        markersize_err_p =  [5] * len(protlist) 
        markerdict = dict(zip(protlist, marker_p))
        sizedict = dict(zip(protlist, markersize_p))
        size_err_dict = dict(zip(protlist, markersize_err_p))
        if "_r" in p:
            prot = p.split("_r")[0]
        else: 
            prot = p
        ax1.errorbar(x, y, yerr = e, fmt='o', markersize=size_err_dict[prot], markeredgecolor='black', linewidth='0.8',capsize=0, c= colordict[prot], label=p, elinewidth = '2')
        ax1.plot(extrapolated_x, extrapolated_y, ls = "--", lw = 0.8,  c= colordict[prot], label='_nolegend_')
        ax2.scatter(x, result.residual, c= colordict[prot], s=sizedict[prot], marker=markerdict[prot], linewidth='0.5', edgecolor='black')

    #set axis boundaries
    ax1.set_xlim([xmin,extended_xmax]) #to see last point
    ax2.set_xlim([xmin,extended_xmax])
    ax1.set_ylim(bottom=0)
    #add minor ticks at every 100 interval 
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(300))
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(300))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax1.tick_params(axis='both') 
    ax2.tick_params(axis='both') 

    #set max number of ticks on y-axis 
    ax2.locator_params(axis='y',nbins=4)

    #outliertext = "Outliers by conc:\n"+str(outliers).strip('[]')+"\n\nInd. Outliers:\n"+str(outliers_ind).strip('[]')
    #ax1.text(extended_xmax*1.05, 1, outliertext, verticalalignment='bottom', rotation="horizontal", )
    fitMatrix = fit.as_matrix()
    handles, labels = ax1.get_legend_handles_labels()
    #lgd = ax1.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1))
    ax1.xaxis.labelpad = 10
    lgd = ax1.legend(handles, labels, loc='upper left', numpoints = 1 )
    fig.show()
    fig.savefig(dirOut + group+'-initialvelocity.pdf', bbox_extra_artists=(lgd,my_suptitle), bbox_inches='tight')

    fit["color_p"] = fit["Protein"].map(colordict)
    fit["label_p"] = fit["Protein"].str.replace('_','\n')
    fout = dirOut + group+'-Kcat_vs_Km.pdf'
    createScatterWithLinearRegression(fit, fout)	

    fit = fit.drop("color_p", axis=1)
    fit = fit.drop("label_p", axis=1)
    fit = fit.set_index("Protein")
    fit = fit.astype('float')

    #temp in kelvin
    T = temp + 273.15 
    #gas constant in kcal/K*mol
    R = 1.985*10**-3
    fit["kcat/km (M-1*s-1)"] = fit["kcat/km (uM-1*s-1)"]*10**-6

    ddG = fit.filter(["kcat (1/s)","kcat_err","kcat/km (uM-1*s-1)","kcat/km_err"], axis=1)
    ddG = ddG.drop('SsWT', axis=0)
    kcat_wt = fit["kcat (1/s)"]["SsWT"] #values 
    ddG["ddg-kcat (kcal/mol)"] = (-R*T*np.log(ddG["kcat (1/s)"]/kcat_wt))
    keff_wt = fit["kcat/km (uM-1*s-1)"]["SsWT"] #values 
    ddG["ddg-Keff (kcal/mol)"] = (-R*T*np.log(ddG["kcat/km (uM-1*s-1)"]/keff_wt))
    ddG2 = ddG.filter(["ddg-kcat (kcal/mol)", "ddg-Keff (kcal/mol)"], axis=1)

    fit = fit.join(ddG2) 
    fit.to_csv(dirOut + group+'-initialvelocity-mm.csv', sep=',' , index=True)



    return 


# In[8]:



extendx = 100

group1 = 'SsMutant-double-pooled'+extracaption
protlist1 = ['SsWT', 'I45A','I45K',  'S70A',  'M73A', 'M73I', 'I45A_S70A', 'I45A_M73A', 'I107K', ]

grouplist = [group1]
protlist = [protlist1]

proteinconc=1
print protlist
for (group, prots) in zip(grouplist, protlist):
    createScatterErrorBarByGroups(velocity, velstd, group, prots, proteinconc)

