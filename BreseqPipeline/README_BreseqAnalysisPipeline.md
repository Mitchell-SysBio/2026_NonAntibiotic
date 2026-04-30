# README_BreseqAnalysisPipeline: 
## The pipeline will take perform breseq pipeline and will use a custom Python script to determine the coverage. Using the breseq results and the coverage information the pipeline will subtract the mutations in the ancestral strain from the condition evolved strains. The pipeline will then compiled the mutations across the replicates evolved clones per condition and map the mutations to genomic coordinates for visualization. The pipline will also determine if the effect on the mutated gene was loss-of-function or modification for mutations that breseq does not call the effect (amplification, new junction, & missing coverage). The effect call is based on the same rules as breseq, which uses the location of the mutation and type of mutation. 


## This pipeline contains files organize into 4 folders:
1.	‘ForRunningBreseq’: will organize sequencing data for running breseq pipeline and will run breseq pipeline in a loop
	- Files:
		- batch_ breseq.sh
		- organizingFastQFiles_V2.m
3.	‘ForOrganizingBreseqOutputs’: will compile and organize breseq output files 
   	- Files:
		- batch_GDrun.sh
		- batch_processData.sh
		- gettingBreseqData.m
		- plotCoverageAllBAM.py 
4.	‘BreseqAnalysisPipeline’: will analyze breseq data and create table outputs and circa plots 
   	- Files:
		- amplificationCompiler_v4.m
		- BreseqAnalyzeWGS_wrapper_V4.m
		- CircaPlotter_v10.m
		- gdColumnNames.mat
		- JCMCcompiler_V2.m
		- predictedMutationCompiler_V6assembly.m
		- Copy-number-and-essentiality-of-all-genes-(E-coli-K-12-substr-MG1655).xlsx

## Requirements
- breseq pipeline installed into a conda environment. Make sure to name the evironment "breseq". See https://github.com/barricklab/breseq/wiki for installation instructions 
- Install prerequisistes for analysis of breseq results and for determining coverage 
    1. Create a conda environment: conda activate breseqAnalysis
    2. Install in environment: 
```
conda install bioconda::pysam
conda install numpy
conda install matplotlib
conda install scipy
pip install html-to-csv
```
- https://github.com/hanwentao/html2csv  

### To run breseq for multiple samples:
1)	In the parent folder put the reference file and the “batch_breseq.sh” script
	    ⁃	Reference file is the reference genome sequence files in GenBank, GFF3, or FASTA format
2)	After moving scripts into parent folder open the metadata file from sequencing company (SeqCoast send us a metadata csv file) and do the following steps.
        Initial steps before running:
        1) Add a column called "fastQFilename" to sample manifest anywhere in the table (the name should be in row 2 which is the header. The first row is just their irrelevant extra information that you won’t extract) (column order doesn’t matter). 
        2) In the column use the concatenate function to combine the “Order ID” and the “SeqCoastTubeID” putting a “_” between them and at the end, and making the SeqCoastTubeID have all the same number of digits eg. 001, 025, 100. For consistency this should always be 3 digits, but it is always best to manually check your samples to make sure this is correct.
            -  To do so type =CONCATENATE(‘Order ID well’, “_”, TEXT(’SeqCoastTubeID’,”000”), “_”)
            -  This will create a name like “OrderNumber_SeqCoastTubeID_” eg 6075_001_
        3) Then drag this formatting to the entire column and save your updated metadata file
            - This is normally the naming convention of the beginning of the SeqCoast FastQ files. You may need to check to make it match the fastq filenames in case they change their naming conventions.
	3)	Use “organizingFastQFiles_V2.m” to move fastQ files into folders and name them as folders(SampleName) 
	    ⁃	This makes a subfolder and puts the R1 and R2 fastq files of one sample into it. The folder is titled with the sample name
	4)	Open the terminal and open the environment with breseq installed: conda activate breseq
	5)	cd to the parent folder 
	6)	Run script in terminal in the parent folder: bash batch_breseq.sh

## To run Breseq Analysis Pipeline:
### Requirements:
•	Parent folder must contain:
        ⁃	Reference file
        ⁃	MATLAB scripts:
            ⁃	ForOrganizingBreseqOutputs/gettingBreseqData.m
Before starting make sure your environment has these two environments installed
	•	Environment with breseq installed: Name = breseq
	•	Environment with html2csv installed: Name = breseqAnalysis

### Step 1: Organize breseq output files
1)	In MATLAB in the parent folder (the directory where all you scripts and subfolders are held) run the script: gettingBreseqData.m
	    ⁃	This will organize the output files into “Results” folder, then 3 sub-folders: htmlfiles, gdfiles, bamfiles
2)	In ‘Results’ folder put bash scripts:
        ⁃	batch_GDrun.sh
        ⁃	batch_processData.sh
3)	In ‘Results’ folder put ‘BreseqAnalysisPipeline’ MATLAB scripts:
        ⁃	amplificationCompiler_v4.m
        ⁃	BreseqAnalyzeWGS_wrapper_V4.m
        ⁃	CircaPlotter_v10.m
        ⁃	gdColumnNames.mat
        ⁃	JCMCcompiler_V2.m
        ⁃	predictedMutationCompiler_V6assembly.m
        ⁃	Copy-number-and-essentiality-of-all-genes-(E-coli-K-12-substr-MG1655).xlsx
4)	In ‘ Results/bamfiles’ folder put ‘plotCoverageAllBAM.py’

### Step 2: Compile summary files from output files
1)	In the terminal, cd to the ‘Results’ folder 
2)	Edit the batch_GDrun.sh: put the name of your reference file 
	    ⁃	Run the script: bash batch_GDrun.sh
	    ⁃	This will activate the ‘breseq’ environment and create a comparison table in ‘html’ and ‘tsv’ format
3)	Run the script: bash batch_processData.sh
	    ⁃	This will activate the ‘breseqAnalysis’ environment 
	    ⁃	Create matrix files with the coverage per base (~30-60 seconds per sample)
	    ⁃	Convert index.html to csv files
4)	Open the ‘AllComparision.html’ in Microsoft Excel
        ⁃	In the columns titled with your strain names:
            1)	Find “~?” and replace with blank
                ⁃	“?” means that coverage was too low to call that mutation, because a small mutation overlaps a missing coverage or junction region
            2)	Copy this delta: “Δ” and find and replace with blank
                ⁃	“Δ” indicates that the mutation in that row is fully contained within a region that is deleted within the specified sample 
            3)	Save as .xlsx file
                ⁃	We are replacing these symbols so that MATLAB can read the comparison table without errors and so the mutation that is fully contained within a region that is deleted within the specified sample isn’t counted as a mutation for the strain 

### Step 3: Analyze breseq data 
- This analysis pipeline will compare across strains evolved on the same drug to the ancestor and create comparison tables and will create circa plots 
	1)	In Results folder put the MAT files from the “BreseqAnalysisPipeline” [insert names], “AllComparision.xlsx”, “ALLGD.tsv” 
	2)	In MATLAB in the ‘Results’ folder, open the ‘BreseqAnalyzeWGS_wrapper_V4.m’
	3)	Fill in the user inputs 
	4)	Run the script 

