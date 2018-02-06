#!/usr/bin/env ipython

# -*- coding: utf-8 -*-



import csv
import pandas as pd
import numpy as np
from numpy.random import normal


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
import itertools
from figure_dict import *


plt.style.use('figures')



dirIn = "input/Function/"
filelist = ["20170628-SsHelix", "20170717-SsHelix-plate1", "20170717-SsHelix-plate2"]

dirOut = "input/Function/cleaned/"



def removeTopRows(df, numRows ):
    df = df.iloc[numRows:]
    df = df.reset_index(drop=True)
    return df

def setRowAsColumnHeader(df, rowNum):
    df.columns = df.iloc[rowNum]
    df = df.reindex(df.index.drop(rowNum))
    return df
   
def trimDataset(df, startrow, endrow):
    numrow = endrow - startrow
    df = df.head(endrow) 
    df = df.tail(numrow) 
    df = df.reset_index(drop=True)
    return df

def get_sample_cond(protein, substrate):	
    samples = []
    for p, s in itertools.product(protein, substrate):
        cond = p+"-"+str(s)+"uMCdRP"
        samples.append(cond)
    return samples

def substrateBackground(df, protein, substrate):
    lst = list(df)
    substrateBackground = pd.DataFrame() 
    for p, s in itertools.product(protein, substrate):
        cond = p+"-"+str(s)+"uMCdRP"
        backSub = df["CdRP-"+str(s)+"uMCdRP"]
        substrateBackground[cond] = backSub
    return substrateBackground

def proteinBackground(df, protein, substrate):
    lst = list(df)
    proteinBackground = pd.DataFrame() 
    for p, s in itertools.product(protein, substrate):
        cond = p+"-"+str(s)+"uMCdRP"
        backProt = df[str(p+"-0uMCdRP")]
        proteinBackground[cond] = backProt
    return proteinBackground

def bufferBackground(df, protein, substrate):
    bufferBackground = pd.DataFrame() 
    for p, s in itertools.product(protein, substrate):
        cond = p+"-"+str(s)+"uMCdRP"
        backBuffer = df["CdRP-0uMCdRP"]
        bufferBackground[cond] = backBuffer
    return bufferBackground
    
def getDF(fin, protein, substrate):
    #dataframe
    df = pd.read_excel(fin)
    df = removeTopRows(df,46)

    #set first row as column header
    df = setRowAsColumnHeader(df, 0)

    #trim dataset (df, startrow, endrow):
    numConditions = (len(protein) * len(substrate)) + 1
    df = trimDataset(df, 1, numConditions)
    df = df.transpose()

    #get sample conditions and set column headers
    sample_cond = get_sample_cond(protein, substrate)
    df.columns = sample_cond 

    #remove temp
    df = df.iloc[1:]
    df.index.name='Time[s]'

    #change data object to float type
    df = df.astype(float)

    #returns dataframes with approprate backgrounds for subtraction 
    subBackground  = substrateBackground(df, protein, substrate)
    protBackground = proteinBackground(df, protein, substrate)
    buffBackground = bufferBackground(df, protein, substrate)

    #subtract background from sample reads 
    df = df.sub(subBackground, axis=0).sub(protBackground, axis=0).add(buffBackground, axis=0)
    df, substrate = removeOutliersByConc(df, substrate, outliers)
    return df

def getVelocity(df, protlist, substrate):
    substrateOnly = "CdRP-"
    df = df.drop([col for col in df.columns if substrateOnly in col],axis=1)
    velocity = pd.DataFrame(columns = ["CdRP[uM]"] + protlist )
    deltaT = df.index.astype(float)
    for j in xrange(len(substrate)):
        vel = []
        vel.append(substrate[j])
        for i in xrange(len(protlist)):
            cond = protlist[i]+"-"+str(substrate[j])+"uMCdRP"
            deltaFl = df[cond]
            #ordinary least square, linear regression  
            slope = sm.OLS(deltaFl,deltaT).fit().params
            vel.append(slope[0])
        velocity.loc[len(velocity)]=vel
    return velocity

def removeOutliersByConc(df, substrate, outliers):
    for o in outliers:
        curr_o = "-"+str(o)+"uM"
        df.drop([col for col in df.columns if curr_o in col],axis=1,inplace=True)
        substrate.remove(o)
    return (df, substrate)

def removeIndividOutliers(velocity, outliers_ind):
    if not outliers_ind:
        print "Individual outliers removed: None\n"
    else: 
        for c in outliers_ind:
            velocity.loc[velocity['CdRP[uM]'] == c[1], c[0]] = np.nan
    return velocity

#-Michalis-Menten 
def min_mm(params, x, y):	
    V = params["V"].value
    K = params["K"].value
    r = params["r"].value
    model = (V * x/(K + x)) - r
    return model - y

def get_mm(x, V, K):	
    Y = (V * x/(K + x))
    return Y


# In[4]:


temp = 30
extracaption = "-"+str(temp)+"C"
numReps = 3
dfs = {}
vels = {}
errorbar = 0


##---variables of conditions---
for f in filelist: 
    substrate = []
    protein = []
    if (f == "20170628-SsHelix"): 
        substrate = [0, 150, 300, 400, 500, 600, 750, 900, 1050, 1200, 1350, 1500]
        protein =   ["SsWT", "I107K", "I45A", "I45A_S70A", "CdRP"]
    elif (f == "20170717-SsHelix-plate1"):
        substrate = [0, 100, 200, 300, 400, 500, 600, 750, 900, 1100, 1300, 1500]
        protein =   ["SsWT", "I45A", "I45K", "S70A", "I45A_S70A", "CdRP"]
    elif (f == "20170717-SsHelix-plate2"):
        substrate = [0, 100, 200, 300, 400, 500, 600, 750, 900, 1100, 1300, 1500]
        protein =   ["M73A", "M73I", "I45A_M73A", "CdRP"]
    else:
        pass
    protlist =  [ p for p in protein if p != "CdRP" ]
    proteinconc = 1 #in units of [uM]
    
    #file names for combined data set
    extracaption = "-combinedreps"
    comb_fout  = dirOut+f+extracaption+'-filtered.csv' 
    comb_fout2 = dirOut+f+extracaption+'-initialvelocity-vi.csv'
    comb_stdev = dirOut+f+extracaption+'-initialvelocity-stdev.csv'


    for i in range(1, numReps+1):
        #samples conditions
        str_i = str(i) #convert int to string 
        #varies with rep 
        outliers = []
        outliers_ind = []
        #files 
        fin = dirIn+f+"-r"+str_i+'.xlsx'  
        fout  = dirOut+f+"-r"+str_i+'-filtered.csv'
        fout2 = dirOut+f+"-r"+str_i+'-initialvelocity-vi.csv'


        #save dataframes
        dfs[i] = getDF(fin, protein, substrate)
        dfs[i].to_csv(fout, sep=',' , index=False)
        vels[i] = getVelocity(dfs[i], protlist, substrate)
        vels[i].to_csv(fout2, sep=',' , index=False)
        #fit = getModel(vels[i], velstd, substrate, protlist)


    dfpanel = pd.Panel(dfs)
    dfmean = dfpanel.mean(axis=0)
    dfstd = dfpanel.std(axis=0)

    velpanel = pd.Panel(vels)
    velmean = velpanel.mean(axis=0)
    vel_stdev = velpanel.std(axis=0)
    velstd = pd.concat([velmean.ix[:,0],vel_stdev.ix[:,1:] ], axis =1)

    dfmean.to_csv(comb_fout, sep=',' , index=False)
    velmean.to_csv(comb_fout2, sep=',' , index=False)
    velstd.to_csv(comb_stdev, sep=',' , index=False)

