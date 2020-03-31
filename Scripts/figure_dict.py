####C:\ProgramData\Anaconda2\Lib\site-packages\matplotlib\mpl-data\stylelib
##print(plt.style.available)


## use \rm to remove italized letters with subscripts
colordictFit = {
'S70A_I107K': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'I45A_S70A': (0.8784313725490196, 0.396078431372549, 0.30980392156862746, 1.0), 
'S70A': (0.8313725490196079, 0.8588235294117647, 0.9019607843137255, 1.0), 
'I45K_I107K': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'I45A_M73A': (0.8313725490196079, 0.8588235294117647, 0.9019607843137255, 1.0), 
'I45K': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'SsWT': (0.8980392156862745, 0.8470588235294118, 0.8196078431372549, 1.0), 
'M73A': (0.4823529411764706, 0.6235294117647059, 0.9764705882352941, 1.0), 
'I45A': (0.8980392156862745, 0.8470588235294118, 0.8196078431372549, 1.0), 
'I45K_I107A': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'I107A': (0.8980392156862745, 0.8470588235294118, 0.8196078431372549, 1.0), 
'I107K': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'I45K_S70A': (0.8, 0.25098039215686274, 0.22745098039215686, 1.0), 
'M73A_I107K': (0.9490196078431372, 0.796078431372549, 0.7176470588235294, 1.0),
'Vector':'grey'}
    
#colordictFit = {
#'S70A_I107K': (0.8661796744841722, 0.23670383185954005, 0.17518902986030963, 1.0), 'I45A_S70A': (0.895317185697808, #0.29799307958477406, 0.20332179930795802, 1.0), 'SsWT': (0.0, 0.4207843137254901, 0.1145098039215691, 1.0), 'I45K_I107K': #(0.7403562732282445, 0.08956555171087942, 0.1508855568371139, 1.0), 'I45A_M73A': (0.0, 0.42823529411764677, #0.04784313725490215, 1.0), 'I45K': (0.8661796744841722, 0.23670383185954005, 0.17518902986030963, 1.0), 'S70A': (0.0, #0.42823529411764677, 0.04784313725490215, 1.0), 'M73A': (0.0, 0.2733333333333334, 0.15111111111111108, 1.0), 'I45A': #(0.06686658977316451, 0.5291708317313859, 0.2788722286300144, 1.0), 'I45K_I107A': (0.7863590926566695, 0.13474817377931464, #0.15264898116109185, 1.0), 'I107A': (0.14053056516724466, 0.6148558246828149, 0.32336793540945824, 1.0), 'I107K': #(0.8322542611815958, 0.1799974368832491, 0.1545508137895681, 1.0), 'I45K_S70A': (0.8322542611815958, 0.1799974368832491, #0.1545508137895681, 1.0), 'M73A_I107K': (0.14053056516724466, 0.6148558246828149, 0.32336793540945824, 1.0), 
#'Vector':'grey'}

colordictCat = {
'SsWT': 'black', 
'S70A_I107K': '#e31a1c', 
'I45A_S70A': '#fdbf6f', 
'S70A':  '#33a02c', 
'I45K_I107K': '#a6cee3', 
'I45A_M73A': '#ffffb3', 
'I45K': '#fb9a99',  
'M73A': '#b15928', 
'I45A': '#6a3d9a', 
'I45K_I107A': '#1f78b4', 
'I107A': '#cab2d6', 
'I107K': '#8dd3c7', 
'I45K_S70A':'#d9d9d9', 
'M73A_I107K': '#ff7f00', 
'Vector':'white', 
}

inch = 2.54
scalefactor=1.5
figsize = {
'Single':8.9/inch*scalefactor,
'Double':18.3/inch*scalefactor,
'Height':24.7/inch*scalefactor
}


#fitColors = [m.to_rgba(x) for x in noVector["SeleCoef"] ]


#'Vector': 'grey',
#'D121A': 'brown', 
#'K207A': 'magenta', 
#'D121A_K207A': 'lavender', 
#'D165A': 'olivedrab',

#'F40A':'coral',
#'D128K':'yellowgreen',
#'V78A':'tan',
#'M73I': 'goldenrod', 
#'L96A':'sienna',
#'A122I':'teal',
#'S70A_A122I':'navy',
#'M73A_A122I':'lightblue',
#'E155A':'royalblue',
#'I107A_E155A':'darkorchid',
#'I107A_K207A':'seagreen',
#'R175A':'lime',
#'D128A_R175A':'darkorange', 
#'D165N':'steelblue', 
#'D165N_R175A':'lightblue', 
#'L108A':'deepskyblue', 
#}

labeldict = {
   
'eff_func_cap':'Effective functional capacity',
'ddG_eff_func_cap':'Relative effective functional capacity',
'RelativeProteinConc':'Normalized [IGPS]',
'selecoef': 's',
'Selcoeff': 's', 
'Selcoeff_all': 's', 
'Selcoef': 's', 
'SeleCoef': 's', 
'doubling_time': 'Doubling time (Hrs)',
'dG_NI': u'$\u0394G_{NI}$'+ u' $ (kcal \cdot mol^{-1})$', 
'dG_IU': u'$\u0394G_{IU}$'+ u' $ (kcal \cdot mol^{-1})$', 
'dG_total': u'$\u0394G_{Total}$'+ r' $ (kcal \cdot mol^{-1})$',
'Km (uM)': u'$K_{m}$' + u'(\u00B5M)',  
'kcat (1/s)': u'$k_{cat}$'+ u' $(s^{-1})$', 
'kcat/km (uM-1*s-1)': u'K$_{eff}$'+ u' $(\u00B5M^{-1} \cdot s^{-1})$',
'kcat/km (nM-1*s-1)': u'K$_{eff}$'+ u' $(nM^{-1} \cdot s^{-1})$', 
'kcat/km (10-3 uM-1*s-1)': u'K$_{eff}$'+ u' $(10^{-3} \cdot \u00B5M^{-1} \cdot s^{-1})$',    
'eff_func_cap': u'$[IGPS]_{i} * K_{eff}$',  
'eff_func_cap_kcat': u'$[IGPS]_{i} * k_{cat}$',
'eff_func_cap_ddGKm': u'$[IGPS]_{i} * \u0394\u0394G_{km}$',   
'eff_func_cap_ddGKcat': u'$[IGPS]_{i} * \u0394\u0394G_{kcat}$',     
'eff_func_cap_ddGkeff': u'$[IGPS]_{i} * \u0394\u0394G_{Keff}$',     
'ProteinExpression': u'normalized [IGPS]'+ u' (R.F.U.)',
'ddG_NI': u'$\u0394\u0394G_{NI}$'+ u' $ (kcal \cdot mol^{-1})$', 
'ddG_IU': u'$\u0394\u0394G_{IU}$'+ u' $ (kcal \cdot mol^{-1})$', 
'ddG_total': u'$\u0394\u0394G_{Total}$',
'ddG_Km (kcal/mol)': u'$\u0394\u0394G_{km}^{\u2021}$', 
'ddG_kcat (kcal/mol)': u'$\u0394\u0394G_{kcat}^{\u2021}$', 
'ddG_keff (kcal/mol)': u'$\u0394\u0394G_{keff}^{\u2021}$', 
'ddG_TS (kcal/mol)': u'$\u0394\u0394G^{\u2021}$'+ u' $(kcal \cdot mol^{-1})$',    
'Tm': u'$T_{m}$'+ u' ($^\circ$C)',
'Temp': u'Temp ($^\circ$C)', 
'RFU-Temp':r'$\Delta RFU \cdot \Delta Temp^{-1}$', 
'MRE': r'$MRE_{222}(10^3 deg \cdot cm^2 \cdot dmol^{-1})$',
'Fapp': 'Fraction unfolded',
'Urea': u'[Urea] (M)',
'Residual': 'Residual',
'interaction_NI': u'$\u03B4_{NI}$' + u'$(kcal \cdot mol^{-1})$',
'interaction_IU': u'$\u03B4_{IU}$' + u'$(kcal \cdot mol^{-1})$',
'interaction': u'\u03B4 ' + u'$(kcal \cdot mol^{-1})$',
'interaction_Km': u'$\u03B4_{km}$',
'interaction_kcat': u'$\u03B4_{kcat}$',
'interaction_keff': u'$\u03B4_{keff}$',
'interaction_s': u'$\u03B4_{s}$',
'wavelength':u'Wavelength (nm)', 
'MRE-spectra':r'$MRE (10^3 deg \cdot cm^2 \cdot dmol^{-1})$', 
'MRE222norm':r'$MRE_{mut}/MRE_{WT}$', 
'CdRP':u'[CdRP] (\u03bcM)', 
'initialvel': u'Initial velocity ($\u0394RFU \cdot s^{-1})$',
'Time':u'Time (Hours)',
'OD600':r'$\rm Log_{2}(OD_{600})$',   
}

axesdict = {
'Selcoeff': [-.2, 0.2, 0.05], 
'SeleCoef': [-.2, 0.2, 0.05], 
'Selcoeff_all': [-.3, 0.3, 0.1], 
'doubling_time': [3, 6, 1],
'dG_NI': [0, 10, 2], 
'dG_IU': [0, 8, 2], 
'dG_total': [6, 18, 2],
'ddG_NI': [-4, 6, 2], 
'ddG_IU': [-4, 1, 1], 
'ddG_total': [-4, 1, 1],
'Km (uM)': [200, 1200, 200], 
'kcat (1/s)': [0, 30, 5], 
'kcat/km (uM-1*s-1)': [0.0001, 0.035, 0.005],
'kcat/km (nM-1*s-1)': [0, 35, 5],
'kcat/km (10-3 uM-1*s-1)': [0, 35, 5],
'ddG_kcat (kcal/mol)': [-0.3, 0.7, 0.2],
'ddG_keff (kcal/mol)': [-0.3, 0.7, 0.2],
'ProteinExpression': [0, 2, 0.5],
'Tm': [45,75, 10],    
}

axesdict2D = {
'Selcoeff': [-.25, 0.15], 
'dG_NI': [1,9], 
'dG_IU': [4,8], 
'dG_total': [7, 15],
'ddG_NI': [-3.5,4.5], 
'ddG_IU': [-2,1], 
'ddG_total': [-4,4],
'Km (uM)': [200, 1200], 
'kcat (1/s)': [0, 30], 
'kcat/km (nM-1*s-1)': [0, 35, 5],
'ProteinExpression': [0, 2, 0.5],
'Tm': [45,70] 
}


orderdict = {
'SsWT': 1, 
'I45A': 2, 
'I45K': 3, 
'S70A': 4,
'M73A': 5, 
'M73I': 6, 
'I107A': 7, 
'I107K': 8,  
'D128A': 9, 
'D121A': 10, 
'D165A': 11,
'K207A': 12, 
'I45A_S70A': 13, 
'I45A_M73A': 14,
'I45K_S70A': 15, 
'I45A_I107A': 16,
'I45K_I107A': 17, 
'I45K_I107K': 18, 
'S70A_I107K': 19,
'M73A_I107K': 20,  
'I107A_D128A': 21, 
'D121A_K207A': 22 
}