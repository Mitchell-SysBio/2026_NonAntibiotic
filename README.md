# 2026_NonAntibiotic

# BreseqAnalysisPipeline
Look at README in that folder for information on breseq analysis pipeline 

# CodeForFigure 
Code to create figures and analyze data

## Code to analyze data for creating figures
### Auto-resistance Folder

### Cross-resistance Folder


### Data files: 
- ABXmechColorMap.mat: contains mechanism and color information for antibiotic classes 
- AllEvo_WGSresults_20260410.mat: contains mutation results of "breseq analysis" 
- allLFC_021926.mat: contains log2 fold change and MGC results of broth microdilution assay 
- AutoResistanceData_2026-04-24.mat: contains IC50 results of dose-response assay 
- FigureOrder.mat: contains order of antibiotic and non-antbiotic names for figures 
- MG1655_NC_00913_gene_positions.csv: contains start and end position of genes in E. coli MG1655 

### Figure 1 & Supplementary Figures 1: 
Figure1.m is a script that outputs Figure 1C, Supplementary Figures 1. It calculates the log2 fold change in IC50 and MGC and plots the bar plots of the log2 fold change.
- Uses data files: 
    - AutoResistanceData_2026-04-24.mat
    - ABXmechColorMap.mat
    - allLFC_021926.mat
    - FigureOrder.mat

### Figure 2 & Supplementary Figures 2: 
Figure2.m is a script that outputs Figure 2 A & B, and Supplementary Figures 2. It plots the cross-resistance and cross-sensitivity to a panel of antibiotics in bar and radar plots. 
- Uses data files:
    - ABXmechColorMap.mat
    - allLFC_021926.mat

### Figure 3 & Supplementary Figure 5: 
Figure3.m is a script that outputs Figure 3 A,B & C and Supplementary Figures 5. It calls on "organizeMutation" function to get the mutated genes for all strains and the position of the mutation. The script plots the type, effect, and position of mutations in pie charts and plots the frequency of mutation of a gene. It also plots the position of the mutated gene in the genome. 
- Uses data files: 
  - ABXmechColorMap.mat
  - AllEvo_WGSresults_20260410.mat
  - FigureOrder.mat 
  - MG1655_NC_00913_gene_positions.csv

### Figure 4: 
Figure4.m is a script that outputs Figures 4 A and B. It calls on "organizeMutation" function to get the mutated genes for all strains and the position of the mutation.The script plots the position of adaptive mutations of the evolved strains on the genome. 
- Uses data files: 
  - ABXmechColorMap.mat
  - AllEvo_WGSresults_20260410.mat
  - MG1655_NC_00913_gene_positions.csv

### Functions
organizeMutation: function will get mutated genes for all strains. It removes some mutations that were found in many strains including the controls and is marginal evidence in the ancestor. It also removes gatZ and lacZ which is deleted in the strain and it's deletion forms a new junction that is captured by breseq. The IS and insertion genes are also removed because they have no set position. 
- Uses data files: 
  - ABXmechColorMap.mat
  - AllEvo_WGSresults_20260410.mat
  - MG1655_NC_00913_gene_positions.csv
