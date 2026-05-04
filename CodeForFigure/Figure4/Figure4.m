%% Figure 4: Adaptive Mutation position on a line 
% Figure 4A: Adaptive Mutation position on a line for all drugs
% Figure 4B: Adaptive Mutation position on a line for multi-drug resistors 
%% Getting data 
cd ../ 
load AllEvo_WGSresults_20260410.mat
load ABXmechColorMap.mat
MG1655Position= readtable("MG1655_NC_00913_gene_positions.csv");

% Organize mutation data 
[allDrugGeneData, superBugs] = organizeMutation(drugGenesTable,MG1655Position, ABXmech, ABX); 
cd Figure4/

uniMech = unique(ABXmech.Mechanism);
rearrangeI = [2,1,8,3,7,4,5,6];
uniMech = uniMech(rearrangeI); % reorganing so it's in desired order 
% Notes:
% allGenes: Some genes have no effect on gene-> shown by NumStrainsLOF = 0 & NumStrainsModification=0
%% Figure 4A: Line plot of mutated genes for all drugs 
% getting drugs for mutated gene at each position: 2>= strains  
allGenes = vertcat(allDrugGeneData{:,1});
allGenes2 = allGenes(allGenes.NumStrainsMutated>=2 & ~(allGenes.NumStrainsLOF==0 & allGenes.NumStrainsModification==0),:);
uniPosition = unique(allGenes2.Start_position);
genesPerPosition = cell(length(uniPosition), 3);
for i1 = 1:length(uniPosition)
    genesPerPosition{i1,1} = uniPosition(i1);
    genesPerPosition{i1,2} = allGenes2.class(allGenes2.Start_position== uniPosition(i1));
    genesPerPosition{i1,3} = 1:length(allGenes2.class(allGenes2.Start_position== uniPosition(i1)));
    genesPerPosition{i1,4} = allGenes2.NumStrainsMutated(allGenes2.Start_position== uniPosition(i1));   
    genesPerPosition{i1,5} = allGenes2.gene(allGenes2.Start_position== uniPosition(i1)); 
end 

fig= figure('Position',get(0, 'Screensize'),'Color','white');
subplot(2,1,1)
hold on 
yline(0, 'k', 'LineWidth',2)
for j1 = 1:height(genesPerPosition)
    colorInd = genesPerPosition{j1,2};
    colorInd(colorInd==0) = 22; % making 0 = index for black 
    mycolors = ABXmech.Colors(colorInd)'; 
    mycolors = cell2mat(mycolors);
    scatter(genesPerPosition{j1,1}, genesPerPosition{j1,3},[], mycolors, 'filled', 'SizeData',100)
    if length(genesPerPosition{j1,3})>1
    text(genesPerPosition{j1,1},-0.2 ,genesPerPosition{j1,5}(1,1),"Rotation",90, 'VerticalAlignment','baseline', 'HorizontalAlignment', 'right')
    end
end 
xlim([0 4641652])
ylim([-1 10])
title('genes with mutations in 2>= strains of drug: colored by drug type')


% Figure 4B: Line plot of mutated genes for multi-drug resistors only 
% getting drugs for mutated gene at each position: 2>= strains  
superbugsInds =zeros(length(superBugs),1);
for sb = 1:length(superBugs)
superbugsInds(sb) = find(contains(string(allDrugGeneData(:,2)), superBugs{sb}));
end 
allGenes = vertcat(allDrugGeneData{superbugsInds,1});
allGenes2 = allGenes(allGenes.NumStrainsMutated>=2 & ~(allGenes.NumStrainsLOF==0 & allGenes.NumStrainsModification==0),:);
uniPosition = unique(allGenes2.Start_position);
genesPerPosition = cell(length(uniPosition), 3);
for i1 = 1:length(uniPosition)
    genesPerPosition{i1,1} = uniPosition(i1);
    genesPerPosition{i1,2} = allGenes2.class(allGenes2.Start_position== uniPosition(i1));
    genesPerPosition{i1,3} = 1:length(allGenes2.class(allGenes2.Start_position== uniPosition(i1)));
    genesPerPosition{i1,4} = allGenes2.NumStrainsMutated(allGenes2.Start_position== uniPosition(i1));   
    genesPerPosition{i1,5} = allGenes2.gene(allGenes2.Start_position== uniPosition(i1)); 
end 

subplot(2,1,2)
hold on;
yline(0, 'k', 'LineWidth',2)
for j1 = 1:height(genesPerPosition)
    colorInd = genesPerPosition{j1,2};
    colorInd(colorInd==0) = 22; % making 0 = index for black 
    mycolors = ABXmech.Colors(colorInd)'; 
    mycolors = cell2mat(mycolors);
    scatter(genesPerPosition{j1,1}, genesPerPosition{j1,3},[], mycolors, 'filled', 'SizeData',100)
    if length(genesPerPosition{j1,3})>=0
        text(genesPerPosition{j1,1},-0.2 ,genesPerPosition{j1,5}(1,1),"Rotation",90, 'VerticalAlignment','baseline', 'HorizontalAlignment', 'right')
    end
end 
xlim([0 4641652])
ylim([-1 10])

title('Superbugs:genes with mutations in 2>= strains of drug: colored by drug type')

% Set the figure's Renderer property to 'painters' before saving as EPS/SVG:
set(fig, 'Renderer', 'painters');
exportgraphics(fig, 'Figure4A_B.eps', 'ContentType', 'vector');

%% Getting number of mutated strains for multi-drug resistors 
multiGenes = allGenes(~(allGenes.NumStrainsLOF==0 & allGenes.NumStrainsModification==0),:); % removing no modification/LOF
uniGenes = unique(multiGenes.gene);
% Count the number of mutated strains for each unique gene
mutatedCounts = zeros(length(uniGenes), 1);
mutatedLOF = zeros(length(uniGenes), 1);
mutatedModification = zeros(length(uniGenes), 1);
uniqueDrugs = zeros(length(uniGenes), 1);
drugs = cell(length(uniGenes), 1);
for k = 1:length(uniGenes)
    mutatedCounts(k) = sum(multiGenes{strcmp(multiGenes.gene, uniGenes(k)), 'NumStrainsMutated'});
    mutatedLOF(k) = sum(multiGenes{strcmp(multiGenes.gene, uniGenes(k)), 'NumStrainsLOF'});
    mutatedModification(k) = sum(multiGenes{strcmp(multiGenes.gene, uniGenes(k)), 'NumStrainsModification'});
    drugs{k} = unique(multiGenes{strcmp(multiGenes.gene, uniGenes(k)), 'drug'});
    uniqueDrugs(k) = length(drugs{k});
end 
mutatedMulti = table(uniGenes,mutatedCounts, mutatedLOF, mutatedModification,uniqueDrugs, drugs);

%% AcrR mutation in multidrug resistors
acrRinfo = allGenes(strcmp(allGenes.gene, 'acrR'),:);
figure
% Pie #2 - position of local mutations - coding / non-coding
codingEffect = sum(acrRinfo.MutationsCoding);
nonCodingEffect = sum(acrRinfo.MutationsIntergenic);
subplot(1,2,1)
piechart([codingEffect, nonCodingEffect], ["coding", "non-coding"])
title(sprintf('%s position of local mutations:\n total mutations = %i\n coding=%i, non-coding=%i', 'acrR', sum([codingEffect, nonCodingEffect]),codingEffect, nonCodingEffect));

% Pie #3 - effect on CDS - loss-of-function / modification
LOFeffect = sum(acrRinfo.MutationsCodingLOF);
modEffect = sum(acrRinfo.MutationsCodingModification);
subplot(1,2,2)
piechart([LOFeffect, modEffect], ["LOF", "modification"])
title(sprintf('%s position of local coding mutations:\n total mutations = %i\n LOF=%i, modification=%i', 'acrR', sum([LOFeffect, modEffect]),LOFeffect, modEffect));

%% lon mutation in multidrug resistors
loninfo = allGenes(strcmp(allGenes.gene, 'lon'),:);