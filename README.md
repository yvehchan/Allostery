# Molecular mechanism of allostery in TIM barrel

## Goal: 
Protein characterization to determine molecular mechanism of allosteric, beneficial mutations found in a trpytophan biosynthetic enzyme called Indole-3-glycerol synthase (IGPS) that folds to the highly prevalent TIM barrel structure. 

## Summary of reseach 
From our primary high throughput saturating mutagenesis fitness screen of 3 orthologous (same protein in different species) TIM-barrel IGPS proteins, a number of mutations were identified to be beneficial, yielding a fitness advantage over the wildtype. Strikingly, some allosteric mutations, distal from the active site, were beneficial. Of particular interest, mutations of  stability elements called beta-alpha hairpin clamps that are found on the opposite face of the protein as the active site, were found to be beneficial. TIM barrel proteins make up about 10% of all known protein structures and these hairpin clamps are conserved in TIM barrel proteins. Integration of analysis from secondary assays to probe structure, stability, and activity will build a working model for the molecular mechanism of these beneficial mutations. This work has implications for protein design and engineering for undustrial purposes.  

## Tools:
  * Python

Scripts for data analysis - fitness, protein expression, secondary structure, stability, activity

## Demo of program

## Main script: master_pipeline.py
  * Calls secondary Python scripts for data analysis

### Structural data analysis 
  * CD-spectra-temp-comparison.py
  
### Stability data analysis 
  * ThermalMelt-spectra
  * TitrationAnalysis-ddg
  
### Fitness data analysis
  * Fitness-avg-scatter
  
### Function/activity data analysis
  * Function-parse-combine-data
  * Function_plot_vi_scatter_keff_vs_kcat

### Relationship between thermal and thermodynamic stability 
  * Tm_vs_dG
  
### Relationship between protein concentration and fitness
  * s_vs_protein_conc
  
### Relationship between protein stability and fitness
  * s_vs_dGNI_Tm
  
### Relationship between protein stability and activity 

### Model for relationship between stabilty, activity, and fitness

