%% Figure 3 & Supplementary Figure 5: Circa Plots
% Figure 3A: Pie Charts comparing position, effect, frequency 
% Figure 3B: circa plots of 2 drugs
% Figure 3C: cira plots of all drugs 
% Supplementary Figure 5: circa plots of controls 
%% Getting data 
cd ../ 
load AllEvo_WGSresults_20260410.mat
load ABXmechColorMap.mat 
load FigureOrder.mat
MG1655Position= readtable("MG1655_NC_00913_gene_positions.csv");

% Organize mutation data 
[allDrugGeneData, superBugs] = organizeMutation(drugGenesTable,MG1655Position, ABXmech, ABX); 
cd Figure3_SuppFig5/
orderDrugs = [ABXorder; nonABXorder];
ABXcompiled = vertcat(allDrugGeneData{[allDrugGeneData{:,3}]==1,1});
nonABXcompiled = vertcat(allDrugGeneData{[allDrugGeneData{:,3}]==2,1});
nonABXcompiled(strcmp(nonABXcompiled.drug,'Streptozotocin'),:) = [];
totCompiled = {ABXcompiled, nonABXcompiled};
%% Figure 3A: Pie Charts comparing 
% Streptozotocin is removed!
% Pie #1 - genomic effect - local/global
% Pie #2 - position of local mutations - coding / non-coding
% Pie #3 - effect on CDS - loss-of-function / modification
f=figure('Position', get(0, 'screensize'));
labels = {'ABX','nonABX'};
r=2; c =3;

for p=1:2
    temp = totCompiled{p};
    % Pie #1 - genomic effect - local/global
    localEffect = sum(temp.MutationsCoding) + sum(temp.MutationsIntergenic); % local effect is only predicted mutation with intergenic/coding calls by breseq- so snps, indels, etc.
    globalEffect = sum(temp.NumMutations)-localEffect;
    subplot(r,c,p)
    piechart([localEffect, globalEffect], ["local", "global"])
    title(sprintf('%s genomic effect:\n total mutations = %i; total strains= %i\n local=%i, global=%i', labels{p},  sum(temp.NumMutations),sum(4*numel(unique(temp.drug))), localEffect, globalEffect))
    
    % Pie #2 - position of local mutations - coding / non-coding
    codingEffect = sum(temp.MutationsCoding);
    nonCodingEffect = sum(temp.MutationsIntergenic);
    subplot(r,c,p+2)
    piechart([codingEffect, nonCodingEffect], ["coding", "non-coding"])
    title(sprintf('%s position of local mutations:\n total mutations = %i\n coding=%i, non-coding=%i', labels{p}, sum([codingEffect, nonCodingEffect]),codingEffect, nonCodingEffect));
    
    % Pie #3 - effect on CDS - loss-of-function / modification
    LOFeffect = sum(temp.MutationsCodingLOF);
    modEffect = sum(temp.MutationsCodingModification);
    subplot(r,c,p+4)
    piechart([LOFeffect, modEffect], ["LOF", "modification"])
    title(sprintf('%s position of local coding mutations:\n total mutations = %i\n LOF=%i, modification=%i', labels{p}, sum([LOFeffect, modEffect]),LOFeffect, modEffect));

end 
exportgraphics(f, 'Figure3A.svg', 'ContentType', 'vector');
%% Figure 3B: Frequency for all genes 
totalCompiled = vertcat(ABXcompiled, nonABXcompiled);
uniGene = unique(totalCompiled.gene); 
frequencyAll = nan(length(uniGene),1); 
for i = 1:length(uniGene)
    frequencyAll(i)= sum(totalCompiled{strcmp(totalCompiled.gene,uniGene{i}), "NumStrainsMutated"});
end 
f1 = figure;hold on;
plot(sort(frequencyAll,'descend'), 'k', 'Linewidth',2)
area(sort(frequencyAll,'descend'))
xlabel('gene')
ylabel('frequency')
xlim([1,length(sort(frequencyAll,'descend'))])
grid on; box on;
title(sprintf('frequency - statistics by gene (%i genes, %i strains)', length(frequencyAll), 40*4))
exportgraphics(f1, 'Figure3B.svg', 'ContentType', 'vector');

%% Figure 3C: Circa plots of 1 antibiotic
% ordered to match radar plots 
spI = 1;
f=figure('Position', get(0, 'screensize'));
cols = gray(4);
cols = cols([4 3 2 1],:);
i = 32;
mutT = allDrugGeneData{i,1};

% 1 circa per strain 
for s1 = 1:4 
    subplot(1,5,s1)
    hold on;
    strainIn = cellfun(@(x) any(x(:)==s1), mutT.StrainMutated); % getting which genes are mutated in which strain 
    % Define the radius and center of the circle
    radius = 100*10; % Change the radius as desired
    center = [0, 0]; % Change the center coordinates as desired
    genomeLength = 4641652;
    % Generate points along the circumference of the circle
    theta = linspace(pi/2,-(3*pi/2), genomeLength); % points around the circle= size of genome
    %theta = linspace(0,2*pi, genomeLength);
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    p = plot(x, y, '-k','LineWidth', 2);
    
    positionMutation = mutT.Start_position;
    positionMutation(positionMutation==0)=1; % since matlab starts with 1
    I1 = find(strainIn == 1);
    scatter(x(positionMutation(I1)), y(positionMutation(I1)),50,'k','filled');
    for I1idx = I1
        text(x(positionMutation(I1idx)), y(positionMutation(I1idx)), mutT{I1idx,"gene"})
    end
    ax = gca;
    axis(ax,'equal');
    axis('off')
    title(s1)
end

% All strains combined 
subplot(1,5,5)
hold on;
% Define the radius and center of the circle
radius = 100*10; % Change the radius as desired
center = [0, 0]; % Change the center coordinates as desired
genomeLength = 4641652;
% Generate points along the circumference of the circle
theta = linspace(pi/2,-(3*pi/2), genomeLength); % points around the circle= size of genome
%theta = linspace(0,2*pi, genomeLength);
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);
p = plot(x, y, '-k','LineWidth', 2);
positionMutation = mutT.Start_position;
positionMutation(positionMutation==0)=1; % since matlab starts with 1
nStrains = mutT.NumStrainsMutated;
[S, I1] = sort(nStrains, 'ascend');
% so white has edge color
scatter(x(positionMutation(I1(S==1))), y(positionMutation(I1(S==1))),50, cols(S(S==1),:),'filled','MarkerEdgeColor',[0 0 0]);
scatter(x(positionMutation(I1(S~=1))), y(positionMutation(I1(S~=1))),50, cols(S(S~=1),:),'filled');
for I1idx = I1'
    text(x(positionMutation(I1idx)), y(positionMutation(I1idx)), mutT{I1idx,"gene"})
end
ax = gca;
axis(ax,'equal');
axis('off')
title(allDrugGeneData{i,2})

exportgraphics(f, 'Figure3C.svg', 'ContentType', 'vector');
%% Figure 3D: Circa plots of all drugs colored by # of strains
% ordered to match radar plots 
drugOrder = [allDrugGeneData(:,2)];
orgInd = nan(length(orderDrugs),1);
for ip = 1:length(orderDrugs)
    temp = strip(orderDrugs{ip},'both',' ');
    temp = erase(temp,' ');
    temp = erase(temp,' ');
    if strcmp(temp, 'Tetracycline HCl')
        orgInd(ip) = find(contains(drugOrder,'Tetracycline'));
    else
        orgInd(ip) = find(contains(drugOrder,temp,'IgnoreCase',true));
    end
end
spI = 1;
f=figure('Position', get(0, 'screensize'));
cols = gray(4);
cols = cols([4 3 2 1],:);
for i = orgInd' 
    mutT = allDrugGeneData{i,1}; 
    subplot(4,10,spI)
    hold on
    spI = spI + 1;
    % Define the radius and center of the circle
    radius = 100*10; % Change the radius as desired
    center = [0, 0]; % Change the center coordinates as desired
    genomeLength = 4641652;
    % Generate points along the circumference of the circle
    theta = linspace(pi/2,-(3*pi/2), genomeLength); % points around the circle= size of genome
    %theta = linspace(0,2*pi, genomeLength);
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    p = plot(x, y, '-k','LineWidth', 2);
    positionMutation = mutT.Start_position;
    positionMutation(positionMutation==0)=1; % since matlab starts with 1
    nStrains = mutT.NumStrainsMutated;
    [S, I1] = sort(nStrains, 'ascend');
    % so white has edge color
    scatter(x(positionMutation(I1(S==1))), y(positionMutation(I1(S==1))),50, cols(S(S==1),:),'filled','MarkerEdgeColor',[0 0 0]);
    scatter(x(positionMutation(I1(S~=1))), y(positionMutation(I1(S~=1))),50, cols(S(S~=1),:),'filled');
    ax = gca;
    axis(ax,'equal');
    axis('off')
    title(allDrugGeneData{i,2})
end
exportgraphics(f, 'Figure3D.svg', 'ContentType', 'vector');

%% Supplementary Figure 5: Circa plots of controls 
% ordered to match radar plots 
spI = 1;
f=figure('Position', get(0, 'screensize'));
cols = gray(4);
cols = cols([4 3 2 1],:);
for i = [17,41] 
    mutT = allDrugGeneData{i,1}; 
    subplot(4,10,spI)
    hold on
    spI = spI + 1;
    % Define the radius and center of the circle
    radius = 100*10; % Change the radius as desired
    center = [0, 0]; % Change the center coordinates as desired
    genomeLength = 4641652;
    % Generate points along the circumference of the circle
    theta = linspace(pi/2,-(3*pi/2), genomeLength); % points around the circle= size of genome
    %theta = linspace(0,2*pi, genomeLength);
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    p = plot(x, y, '-k','LineWidth', 2);
    positionMutation = mutT.Start_position;
    positionMutation(positionMutation==0)=1; % since matlab starts with 1
    nStrains = mutT.NumStrainsMutated;
    nStrains(nStrains==5) = 4; % making 5 strains into 4
    [S, I1] = sort(nStrains, 'ascend');
    % so white has edge color
    scatter(x(positionMutation(I1(S==1))), y(positionMutation(I1(S==1))),50, cols(S(S==1),:),'filled','MarkerEdgeColor',[0 0 0]);
    scatter(x(positionMutation(I1(S~=1))), y(positionMutation(I1(S~=1))),50, cols(S(S~=1),:),'filled');
    for I1idx = I1'
        text(x(positionMutation(I1idx)), y(positionMutation(I1idx)), mutT{I1idx,"gene"})
    end
    ax = gca;
    axis(ax,'equal');
    axis('off')
    title(allDrugGeneData{i,2})
end

exportgraphics(f, 'Figure3supplementalControl.svg', 'ContentType', 'vector');






