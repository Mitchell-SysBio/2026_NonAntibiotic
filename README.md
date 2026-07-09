# 2026_NonAntibiotic

# BreseqAnalysisPipeline Folder: 
Breseq analysis pipeline will run breseq (https://github.com/barricklab/breseq), determine amplifications, and compile together mutations for replicates of each drug condition. More details are in the README in the folder. The output is "AllEvo_WGSresults_20260410.mat".

# CodeForFigure Folder: 
Code to create figures and analyze data

# Code to analyze data for creating figures
## Auto-resistance Folder: all code and data to determine auto-resistance
### Code:
- gettingDoseResponseData.m: script to organize raw dose response data from "GrowthCurve_All_Results.xlsx" into "Evo_IC50Data_2026-04-22.mat" using the metadata from "AllGrowthCurveMetadata_V2.xlsx".
- AutoResistance_ICdata.m: script to calculate the IC50 for each drug and strain at user-inputed timepoint. Uses the dose response data from "Evo_IC50Data_2026-04-22.mat". Outputs the IC50 into "AutoResistanceData_2026-04-24.mat", and graphs of normalized growth curves and the dose-response curves. Uses the "calcDoseResponse_v2.m" function to calculate the IC50.
- calcDoseResponse_v2.m: function to calculate the IC50 at a user specified time point
### Data: 
- GrowthCurve_All_Results.xlsx: each sheet corresponds to one plate of a Logphase dose-response experiment
- AllGrowthCurveMetadata_V2.xlsx: contains data on which sheet and cell in "GrowthCurve_All_Results.xlsx" corresponds to a drug, drug concentration, and strain
- Evo_IC50Data_2026-04-22.mat: raw dose-response data now organized for each drug, concentration, and strain

## Cross-resistance Folder: all code and data to determine cross-resistance 
### Code
- organizeMGCdata_FINAL.m: script to organize the MGC results and calculate the log2 fold change from the ancester. The MGC results were determined using the AssiST Pipeline (https://github.com/Mitchell-SysBio/AssiST) and are in mat files. The script takes the mat files from the different sets of broth microdilution experiments and correlates it with the strain and drug information from the metadata files ('platemap..xlsx'). It then determines the median MGC of all the ancestor replicates for each antibiotic and uses that determine the log2FC. The script will output the log2FC in "allLFC_021926.mat". 
### Data: 
- Mat files of MGC results: MGCresults4evo_05-Sep-2024.mat, MGCresults5FU_RFB_26-Sep-2024.mat, MGCresultsABX_06-Sep-2024.mat, MGCresultsCarmofur_02-Jul-2025.mat,MGCresultsChlorhexidine_12-Sep-2024.mat, MGCresultsNonABX_06-Sep-2024.mat
- Metadata files: platemap4Evo.xlsx, platemapABX.xlsx, platemapCarmofurRepeat.xlsx, platemapNonABX.xlsx

# Code and data for creating figures 
### Data files: 
- ABXmechColorMap.mat: contains mechanism and color information for antibiotic classes 
- AllEvo_WGSresults_20260410.mat: contains mutation results of "breseq analysis pipeline" 
- allLFC_021926.mat: contains log2 fold change and MGC results of broth microdilution assay 
- AutoResistanceData_2026-04-24.mat: contains IC50 results of dose-response assay 
- FigureOrder.mat: contains order of antibiotic and non-antbiotic names for figures 
- MG1655_NC_00913_gene_positions.csv: contains start and end position of genes in E. coli MG1655 

### Figure 1 & Supplementary Figures 1 & 2: 
Figure1.m is a script that outputs Figure 1C, Supplementary Figures 1 & 2. It calculates the log2 fold change in IC50 and MGC and plots the bar plots of the log2 fold change. It compares the IC50 and MGC using box plots.
- Uses data files: 
    - AutoResistanceData_2026-04-24.mat
    - ABXmechColorMap.mat
    - allLFC_021926.mat
    - FigureOrder.mat

### Figure 2 & Supplementary Figures 3: 
Figure2.m is a script that outputs Figure 2 A & B, and Supplementary Figures 2. It plots the cross-resistance and cross-sensitivity to a panel of antibiotics in bar and radar plots. 
- Uses data files:
    - ABXmechColorMap.mat
    - allLFC_021926.mat

### Figure 3: 
Figure3.m is a script that outputs Figure 3 A,B & C. It calls on "organizeMutation" function to get the mutated genes for all strains and the position of the mutation. The script plots the type, effect, and position of mutations in pie charts and plots the frequency of mutation of a gene. It also plots the position of the mutated gene in the genome. 
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
