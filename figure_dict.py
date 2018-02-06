
colordict = {
'SsWT': 'black', 
'I45A': 'red', 
'I45K': 'tomato', 
'S70A': 'fuchsia',
'M73A': 'yellow', 
'M73I': 'goldenrod', 
'I45A_S70A': 'green', 
'I45A_M73A': 'blue',
'I107A': 'purple', 
'I45A_I107A': 'orange',
'D128A': 'cyan', 
'I107A_D128A': 'darkseagreen', 
'D121A': 'brown', 
'K207A': 'magenta', 
'D121A_K207A': 'lavender', 
'D165A': 'olivedrab',
'I107K': 'indigo',
'F40A':'coral',
'D128K':'yellowgreen',
'V78A':'tan',
'L96A':'sienna',
'A122I':'teal',
'S70A_A122I':'navy',
'M73A_A122I':'lightblue',
'E155A':'royalblue',
'I107A_E155A':'darkorchid',
'I107A_K207A':'seagreen',
'R175A':'lime',
'D128A_R175A':'darkorange', 
'D165N':'steelblue', 
'D165N_R175A':'lightblue', 
'L108A':'deepskyblue' 
}

labeldict = {
'Selcoeff': 's', 
'Selcoeff_all': 's', 
'doubling_time': 'Doubling time (Hrs)',
'dG_NI': u'$\u0394G_{NI}$'+ u' $ (kcal \cdot mol^{-1})$', 
'dG_IU': u'$\u0394G_{IU}$'+ u' $ (kcal \cdot mol^{-1})$', 
'dG_total': u'$\u0394G_{Total}$'+ r' $ (kcal \cdot mol^{-1})$',
'kcat (1/s)': u'$k_{cat}$'+ u' $(s^{-1})$', 
'kcat/km (uM-1*s-1)': u'K$_{eff}$'+ u' $(\u00B5M^{-1} \cdot s^{-1})$',
'kcat/km (10-3 uM-1*s-1)': u'K$_{eff}$'+ u' $(10^{-3} \cdot \u00B5M^{-1} \cdot s^{-1})$',       
'ProteinExpression': u'normalized [IGPS]'+ u' (R.F.U.)',
'ddG_NI': u'$\u0394\u0394G_{NI}$'+ u' $ (kcal \cdot mol^{-1})$', 
'ddG_IU': u'$\u0394\u0394G_{IU}$'+ u' $ (kcal \cdot mol^{-1})$', 
'ddG_total': u'$\u0394\u0394G_{Total}$',
'ddg-kcat (kcal/mol)': u'$\u0394\u0394G_{kcat}^{\u2021} (s^{-1})$', 
'ddg-Keff (kcal/mol)': u'$\u0394\u0394G_{Keff}^{\u2021} (s^{-1})$', 
'Tm': u'$T_{m}$'+ u' ($^\circ$C)',
'Temp': u'Temp ($^\circ$C)', 
'RFU-Temp':r'$\Delta RFU \cdot \Delta Temp^{-1}$', 
'MRE': r'$MRE_{222}(10^3 deg \cdot cm^2 \cdot dmol^{-1})$',
'Urea': u'[Urea] (M)',
'Residual': 'Residual',
'interaction_NI': u'$\u03B4_{NI}$' + u'$(cal \cdot mol^{-1})$',
'interaction_IU': u'$\u03B4_{IU}$' + u'$(cal \cdot mol^{-1})$',
'interaction': '\u03B4' + u'$(cal \cdot mol^{-1})$',
'wavelength':u'Wavelength (nm)', 
'MRE-spectra':r'$MRE (10^3 deg \cdot cm^2 \cdot dmol^{-1})$',   
'CdRP':u'[CdRP] (\u03bcM)', 
'initialvel': u'Initial velocity (\u0394RFU/s)'
}

axesdict = {
'Selcoeff': [-.15, 0.15, 0.05], 
'Selcoeff_all': [-.3, 0.3, 0.1], 
'doubling_time': [3, 6, 1],
'dG_NI': [0, 10, 2], 
'dG_total': [6, 18, 2],
'kcat (1/s)': [0, 30, 5], 
'kcat/km (uM-1*s-1)': [0.0001, 0.035, 0.005],
'kcat/km (10-3 uM-1*s-1)': [0, 35, 5],
'ProteinExpression': [0, 2, 0.5],
'Tm': [45,75, 10] 
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
'I45A_I107A': 15,
'I107A_D128A': 16, 
'D121A_K207A': 17 
}
