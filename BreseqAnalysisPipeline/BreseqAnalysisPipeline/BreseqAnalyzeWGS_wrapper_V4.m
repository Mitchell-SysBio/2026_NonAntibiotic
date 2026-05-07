%% Analyze WGS data from Breseq Results - Wrapper script 
% Updated 2025-01-08 by Carmen Li

% USER INPUTS:
% input name of COMPARE Excel file: 
% NEED TO ADD IN MORE DETAILED INSTRUCTIONS 
% Must replace ?=50%, Δ= 25% for strain
% results before saving as ".xlsx" (these will be treated as a mutation in
% the gene in the circa plots) 
compareFile = "AllComparison.xlsx"; 

% input name of ancestor file: must have an ancestor file
ancestorName = "MG1655_lacZ23_A1_Ancestor_TubeID19"; 

% Results will be compared: Sample name will be split by "_" & compared
% using part of sample name
% EXAMPLE: "A10_Doxifluridine_C1_TubeID97" is compared using 'Doxifluridine' 
% "A10_Doxifluridine_C1_TubeID97" split = ["A10", "Doxifluridine", "C1",
% "TubeID97"] so splitNumber = 2
splitNumber = 2;

% input length of genome 
genomeLength = 4641652; 

% Output mat file name:
outputName = "AllEvo_WGSresults_20260410.mat";

% input name of GD TSV file: made using gdtools COMPARE 
gdFile = "ALLGD.tsv";

% BAM folder name: contains '.mat' file of each strain with the name being
% 'strainname.mat' 
BAMfolder = 'bamfiles';

% HTML folder name: contains 'index.html' of each strain with the name
% being 'strainname.html' 
htmlFolder = 'htmlfiles'; 

% Write over saved files if they already exist?
% 0 = do not write over; 
% 1 = write over 
% 2= plots/ excels aren't created to save memory/time
% 3= only save Circa information plots (saves a lot of time) 
noWriteOver = 3; 

% Amplification parameters: User can change these, but these values seem to be good 
locS = 5*10^4; % dimension of local smoothing
gloS = 2*10^6; % dimension of global smoothing
covCutoff=0.4; % fold cutoff for determining large amplification
ampLengthCutoff = 86; % cutoff for determining length of amplication: anything shorter does not count, right now based on smallest gene in E. coli Kpdf (29 amino acids) 
smallCovCutoff = 5; % fold cutoff for determining small amplification
%% Making directory for saving 

if exist('./CircaPlots', 'dir') ~= 7
    mkdir('CircaPlots')
end

% making directory for Circa plots 
if exist('./CircaInformation', 'dir') ~= 7
    mkdir('CircaInformation')
end

% making directory for Circa plots as figure files
if exist('./CircaFigure', 'dir') ~= 7
    mkdir('CircaFigure')
end

% making directory for excel files
if exist('./BreseqInformation', 'dir') ~= 7
    mkdir('BreseqInformation')
end

% making directory for Amplification plots 
if exist('./AmplificationPlots', 'dir') ~= 7
    mkdir('AmplificationPlots')
    mkdir('AmplificationPlots/NoAmplification')
end

%% Compiling data fom source files 
% ~ Predicted mutations: Getting data from COMPARE & TSV file and remove ancestor mutation 
fprintf('Getting predicted mutations data\n')
load("gdColumnNames.mat")
[compareData] = predictedMutationCompiler_V6assembly(compareFile, gdFile, ancestorName, gdColumnNames);
fprintf('Finished getting predicted mutation data\n')

% ~ Getting amplification data: 
% coverage per base determine using 'plotCoverageAllBAM.py' 
fprintf('Running coverage calculations \n')
[amplificationData] = amplificationCompiler_v4(locS, gloS, covCutoff, BAMfolder, ancestorName, ampLengthCutoff, smallCovCutoff);
fprintf('Finished getting amplifications\n')
cd ./ % returning to main folder

% ~ Importing in junction & missing coverage evidence and saving for each drug
% Using CSV files from html2csv

% Getting files 
fprintf('Getting JC and MC data \n')
[JC, MC, ancJC,ancMC] = JCMCcompiler_V2(ancestorName, htmlFolder);
fprintf('Finished getting JC and MC data \n')

save(outputName) % saving all variables

%% Plotting Circa plot where the mutations are for each strain for each drug 
% getting indexes of strains evolved on same drug
fprintf('Plotting Circa Plots \n')

variables = compareData.Properties.VariableNames; % This was added because reference assemblies have an additional seq_id column at the beginning in compare table
mutIndex = find(strcmp(variables, 'mutation')); % Finding the column that has mutation in order to know how to index 
mutIndexA = find(strcmp(variables, 'annotation')); % Finding the column that has annotation in order to know how to index 
colNames = compareData.Properties.VariableNames(mutIndex+1:mutIndexA-1);
colNames(strcmp(colNames, ancestorName)) = []; % removing ancestor
strainNames = cellfun(@(x) strsplit(x,"_"), colNames, 'UniformOutput',false);
colSplit = cellfun(@(x) x(splitNumber),strainNames);
typeCMP = unique(colSplit); % conditions to be compared
Ssize = get(0, 'Screensize'); % getting screensize
radMult = 100; % radius multiplier to change size of circles 
drugGenes = cell(length(typeCMP),1); % holds information for strains evolved in condition 

for d = 1:length(typeCMP) 
    drugData = CircaPlotter_v10(compareData, genomeLength, JC, MC, Ssize,radMult, amplificationData, typeCMP{d}, noWriteOver, colSplit, splitNumber);
    drugGenes{d} = drugData;
end
temp = cell2table(drugGenes);
drugGenesTable = splitvars(temp,"drugGenes","NewVariableNames", ["CommonGene","PredictedMutationTable","JunctionTable", "MissingCoverageTable", "AmplificationTable", "AllMutatedGeneTable", "StrainInfo"]);
drugGenesTable = addvars(drugGenesTable, typeCMP', 'Before','CommonGene','NewVariableNames','DrugName');

save(outputName, 'drugGenes', 'drugGenesTable', '-append')
fprintf('Finished!\n')