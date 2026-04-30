function [drugGenes] = CircaPlotter_v10(compareData, genomeLength, JC, MC, Ssize,radMult, amplificationData, drugName, noWriteOver, colSplit, splitNumber)
%% Updated 2026-04-01 by Carmen Li: corrected error where number of
% Predicted strains was incorrect because did not subtract AMP from type 
% Plots circa plots and compiles all information for the condition
% Replacing with empty strings for gene name columns that don't have anything
% Added in counting # of mutations per drug for the different categories
% because it's possible to have multiple mutations in the same gene in the
% same strain 
if strcmp(class(compareData.genes_inactivated), 'double')
    compareData.genes_inactivated = repmat("",height(compareData),1);
end
if strcmp(class(compareData.genes_overlapping), 'double')
    compareData.genes_overlapping = repmat("",height(compareData),1);
end
if strcmp(class(compareData.genes_promoter), 'double')
    compareData.genes_promoter = repmat("",height(compareData),1);
end

if noWriteOver~=2 % not skipping plotting
    % creating uifigure
    f = uifigure('color','w','Position',[1 1 Ssize(3) Ssize(4)], 'Visible', 'off');
    a = uiaxes(f, 'Position', [0 0 Ssize(3)*0.4 Ssize(4)],'Visible','off');
    hold(a, 'on');
    title(a, drugName, 'FontSize',16, 'Visible', 'on')
end
drugInds1 = find(strcmp(colSplit,drugName));
variables = compareData.Properties.VariableNames; % This was added because reference assemblies have an additional seq_id column at the beginning in compare table
mutIndex = find(strcmp(variables, 'mutation')); % Finding the column that has mutation in order to know how to index

colInds = drugInds1 + mutIndex; % getting column numbers of samples

%drugInds = find(contains([JC(:,splitNumber+1)],drugName)); % getting index for JC
try
drugInds = find(contains(JC(:,splitNumber+1),drugName)); % getting index for JC
catch
    drugInds = find(strcmp([JC{:,splitNumber+1}],drugName)); % getting index for JC
end 

% intitialzing variables
colors = lines(length(colInds));
legplots = length(colInds);
ampPositions = cell(length(colInds),1);
AstrainInd = cell(length(colInds),1);
drugAMP = cell(length(colInds),1);
mutPositions = cell(length(colInds),1);
mutStrain = cell(length(colInds),1);
mutType = cell(length(colInds),1);
mutClass = cell(length(colInds),1);
mutGenes =  cell(length(colInds),1);
mutFullName = cell(length(colInds),1);
mutFullInfo = cell(length(colInds),1);
drugJC = cell(length(colInds),1);
drugMC = cell(length(colInds),1);


% looping through strains for the condition
for j = 1:length(colInds) % looping through samples
    % Define the radius and center of the circle
    radius = radMult+j*10; % Change the radius as desired
    center = [0, 0]; % Change the center coordinates as desired

    % Generate points along the circumference of the circle
    thetaV = linspace(pi/2,-(3*pi/2), genomeLength); % points around the circle= size of genome
    %theta = linspace(0,2*pi, genomeLength);
    x = center(1) + radius * cos(thetaV);
    y = center(2) + radius * sin(thetaV);

    if j == 1 % storing first circle
        firstCircle = [x',y'];
    end
    if noWriteOver~=2 % not skipping plotting
        % Plot the circle
        p = plot(a, x, y, 'LineWidth', 2, 'Color',colors(j,:));
        legplots(j) = p;
    end
    % getting amplification data
    sName = compareData.Properties.VariableNames(colInds(j));
    aInd = find(contains(amplificationData(:,1), sName));

    if ~isempty(aInd)
        AmpData = amplificationData{aInd,2}; % all positions of amplification
        AmpDataPos = amplificationData{aInd,3}; % positions of amplification seperated by each amplification
        strainAMP = amplificationData(aInd,:);
        if noWriteOver~=2 % not skipping plotting
            % Plot thicker line of amplification
            for ad = 1:height(AmpDataPos)
                pos = AmpDataPos(ad,1): AmpDataPos(ad,2);
                if length(pos)<10000
                    % Plotting small amplifications
                    plot(a,x(AmpDataPos(ad,1): AmpDataPos(ad,2)),y((AmpDataPos(ad,1): AmpDataPos(ad,2))),"-", 'LineWidth', 20, 'Color',colors(j,:));
                else
                    % plotting large amplifications
                    plot(a,x(AmpDataPos(ad,1): AmpDataPos(ad,2)),y((AmpDataPos(ad,1): AmpDataPos(ad,2))), 'LineWidth', 5, 'Color',colors(j,:));
                end
            end
            ampPositions{j} = AmpData;
        end 
     else
            strainAMP = [];
    end
    % Plot mutated genes on the circle
    strainInds = find(~isnan(compareData{:,colInds(j)})); % getting genes that strains mutated
    % assigning color based on mutation class
    fillColor= zeros(length(strainInds), 3)+150/255; % making grey the default = modification
    fillColor(strcmp(compareData{strainInds,"MutationClass"}, "LOF"),:) = 0; % black
    fillColor(strcmp(compareData{strainInds,"MutationClass"}, "noImpact"),:) = 1; % white
    try
        fillColor(strcmp(compareData{strainInds,"MutationClass"}, "AMP"),:) = colors(j,:); % matching line
    catch
        % do nothing
    end
    positionMutation = compareData{strainInds, "position"}; % getting position of the gene
    AstrainInd{j} = strainInds; % storing strain inds
    if noWriteOver~=2 % not skipping plotting
        scatter(a,x(positionMutation), y(positionMutation),100, fillColor,'filled','MarkerEdgeColor',[150/255 150/255 150/255]);
    end

    strainName = compareData.Properties.VariableNames(colInds(j));
    %% Getting new junction evidence
    strainJC = JC{drugInds(j),4};
    geneJC = []; % individual name of genes
    FullGeneJC = [];% fullname of the gene: geneA_geneB
    posGeneJC = [];
    JC_StrainID =cell(height(strainJC),1); %strain name
    mutClassJC = []; % Class: LOF, modification
    for jc = 1:height(strainJC)
        if ~contains(strainJC{jc, "gene_A"},'IS') || ~contains(strainJC{jc, "gene_A"},'ins') % checking if it is NOT an IS element
            posJC = str2double(regexp(strainJC{jc, "position_A"},'\d*','Match'));
            fullGname = sprintf("%s_%s",strainJC{jc,"gene_A"},strainJC{jc,"gene_B"});
        else
            posJC = str2double(regexp(strainJC{jc, "position_B"},'\d*','Match'));
            fullGname = sprintf("%s_%s",strainJC{jc,"gene_B"},strainJC{jc,"gene_A"});
        end
        genesss = [strainJC{jc, 'LOFgenes'}{1}, strainJC{jc, 'ModifiedGenes'}{1}];
        geneJC= [geneJC, genesss];
        FullGeneJC = [FullGeneJC; repmat({fullGname}, length(genesss), 1)];
        mutClassJC = [mutClassJC; repmat("LOF",length(strainJC{jc, 'LOFgenes'}{1}),1); repmat("modification", length(strainJC{jc, 'ModifiedGenes'}{1}),1)];
        posGeneJC = [posGeneJC; zeros(length(genesss),1)+posJC];
        if noWriteOver~=2 % not skipping plotting
            if ~isempty(length(strainJC{jc, 'LOFgenes'}{1}))
                % Has 1 gene with LOF
                plot(a,x(posJC), y(posJC), '*', 'MarkerSize',10, 'LineWidth',2, 'Color', 'k')
            else
                % Modification
                plot(a,x(posJC), y(posJC), '*', 'MarkerSize',10, 'LineWidth',2, 'Color', zeros(1,3)+150/255)
            end
        end
        JC_StrainID(jc) = strainName;
    end
    strainJC = addvars(strainJC,JC_StrainID,'Before','SeqID_A'); % adding strain ID
    if height(strainJC)~=0
        mutTypeJC = cell(length(geneJC),1);
        mutTypeJC(:) = {'JC'};
    else
        mutTypeJC = {};
    end

    %% Getting missing coverage
    strainMC = MC{drugInds(j),4};
    geneMC = {}; lp=1;
    posGeneMC = [];fullNameMC = cell(1000,1);
    MC_StrainID =cell(height(strainMC),1); %strain name
    mutClassMC = {};
  
    for mc = 1:height(strainMC)
        % handling if the side position is multiple locations ex. 1-2
        mcS = str2double(strainMC{mc, 'start'});
        mcE = str2double(strainMC{mc, 'end'});
        if isnan(mcS)
            numS = strsplit(strainMC{mc, 'start'},'–');
            mcS = int64(mean(str2double(numS{1}), str2double(numS{2}))); % getting mean of the positions and making so its an integer
        end
        if isnan(mcE)
            numS = strsplit(strainMC{mc, 'end'},'–');
            mcE = int64(mean(str2double(numS{1}), str2double(numS{2}))); % getting mean of the positions and making so its an integer
        end
        mcPos = mcS:mcE;
        if noWriteOver~=2 % not skipping plotting
            plot(a,x(mcPos) , y(mcPos), 'square', 'MarkerSize', 10,'LineWidth',2,'Color', 'k');
        end
        % getting all genes covered by deletion
        if sum(strfind(strainMC{mc,"gene"}{1,1}, '/')) >0 % these are in the intergenic region
            manyGenes = strainMC{mc,"gene"}{1,1};
            tgene = strsplit(manyGenes,{'/'});
            mutClassMC = [mutClassMC; repmat("modification",length(tgene),1)];
        elseif contains(strainMC{mc,"gene"}{1,1},'–') % many genes 
            manyGenes = strainMC{mc,"description"}{1,1};
            manyGenes(1:(4+strfind(manyGenes,'genes')))=[]; % removing prior to where description says genes 
            tgene = strsplit(manyGenes,{','});
            mutClassMC = [mutClassMC; repmat("LOF",length(tgene),1)];
        else % If only 1 gene 
            tgene = strainMC{mc,"gene"};
            mutClassMC = [mutClassMC; repmat("LOF",length(tgene),1)];
        end 
        tgene = strtrim(erase(tgene, ["[", "]", "→", "←",char(160),' '])); % stripping of special characters
        geneMC = [geneMC,tgene];
        posGeneMC = [posGeneMC; zeros(length(tgene),1)+double(round(median(mcPos)))]; %getting middle of missing coverage
        MC_StrainID(mc) = strainName;
        fullNameMC(lp:lp-1+length(tgene))= {erase(strainMC{mc,"gene"}, [ "→", "←",char(160),' '])};
        lp = lp+length(tgene);
    end
    strainMC = addvars(strainMC,MC_StrainID,'Before','seq_id'); % adding strain ID
    fullNameMC(length(posGeneMC)+1:end) = [];
    mutTypeMC = {};
    if height(strainMC)~=0
        mutTypeMC = cell(length(posGeneMC),1);
        mutTypeMC(:) = {'MC'};
    end

    %% Getting predicted mutations and handling large deletion
    
    genePM = {}; posGenePM = []; mutTypePM = []; fullNamePM=[];mutClassPM = {}; fullInfoPM = []; mutPositionTypePM = {}; 
    for pm = strainInds'
        LOF = {}; MOD = {}; AMP = {}; NIG = {};
        if ~isempty(compareData{pm,'LOFgenes'}{1})
            LOF = compareData{pm,'LOFgenes'}{1};
            if isempty(LOF{1})
                LOF = '';
            end 
        end
        if ~isempty(compareData{pm,'Modifiedgenes'}{1})
            MOD = compareData{pm,'Modifiedgenes'}{1};
            if isempty(MOD{1})
                MOD = '';
            end 
        end
        if ~isempty(compareData{pm,'AmplifiedGenes'}{1})
            AMP = strsplit(compareData{pm,'AmplifiedGenes'}{1},',');
        end
        if ~isempty(compareData{pm,'noImpactGenes'}{1})
            NIG = strsplit(compareData{pm,'noImpactGenes'}{1},',');
        end
        tgene = [LOF, MOD, AMP, NIG];
        tgene = tgene(~cellfun('isempty', tgene)); % removing empty cells
        
        % Getting coding or intergenic for mutation
        if isempty(compareData{pm,"gene_position"}{1,1})
            mutPosType = "NA"; % mutations that are not truly local 
        elseif contains(compareData{pm,'gene_position'}{1,1}, 'intergenic')
            mutPosType = "intergenic";
        else 
            mutPosType = "coding";
        end 
        genePM= [genePM, tgene];
        posGenePM = [posGenePM; zeros(length(tgene),1)+compareData{pm,'position'}];
        mutClassPM = [mutClassPM; repmat("LOF", length(LOF), 1); repmat("modification", length(MOD), 1); repmat("amplification", length(AMP), 1); repmat("noImpact", length(NIG), 1)];
        mutTypePM = [mutTypePM; repmat(mutPosType, length(tgene),1)];%[mutTypePM; repmat({compareData{pm,'mutation_category'}}, length(tgene),1)];
        fullNamePM = [fullNamePM; repmat({compareData{pm,'gene'}{1}}, length(tgene),1)];
        fullInfoPM = [fullInfoPM; repmat({compareData{pm,'mutation'}{1}}, length(tgene),1), repmat({compareData{pm,'annotation'}{1}}, length(tgene),1),repmat({compareData{pm,'description'}{1}}, length(tgene),1)];
    end

    % Saving mutated genes and positions
    mutGenes{j} = horzcat(genePM, geneJC, geneMC);
    posInfo = vertcat(posGenePM, posGeneJC, posGeneMC);
    mutPositions{j} = posInfo;
    mutStrain{j} = zeros(length(posInfo),1) + j;% getting which strain corresponds to mutated gene
    mutType{j} = vertcat(mutTypePM,mutTypeJC,mutTypeMC); % getting type of mutation
    mutClass{j} = vertcat(mutClassPM, mutClassJC, mutClassMC); % getting class of mutation: LOF, modification, amplification
    mutFullName{j} = vertcat(fullNamePM, FullGeneJC, fullNameMC);
    mutFullInfo{j} = vertcat(fullInfoPM, cell(length(FullGeneJC),3), cell(length(fullNameMC),3)); % getting full information 
    drugJC{j} = strainJC;
    drugMC{j} = strainMC;
    drugAMP{j} = strainAMP;
    if noWriteOver~=2 % not skipping plotting
        axis(a,'equal'); % Make the aspect ratio equal to see a perfect circle
    end
end
%% converting into more accessible format
mutPositions = cell2mat(vertcat(mutPositions));
mutStrain = cell2mat(mutStrain);
mutGenes = string(horzcat(mutGenes{:}))';
mutType = string(vertcat(mutType{:}));
mutClass = string(vertcat(mutClass{:}));
mutFullName = string(vertcat(mutFullName{:}));
mutFullInfo = vertcat(mutFullInfo{:});
drugJC = vertcat(drugJC{cellfun(@(x) ~isempty(x), drugJC)});
drugMC = vertcat(drugMC{cellfun(@(x) ~isempty(x), drugMC)});

%% Amplifications: plotting radian lines of mutated genes that are in amplification
% Find genes that are in amplification region by looking at all genes
% mutated for all strains
% We are interested in those that overlap with amplification region
% 2026-04-01: Edited it so only getting unique(mutGenes(loc1)) so that even
% if multiple other strains have a mutated gene that overlaps with the
% amplification the mutated gene is only counted 1 time
for o2 = 1:height(ampPositions) % in order of strain from 1 to 4
    if ~isnan(ampPositions{o2}) %only strains with amplication
        loc1 = find(ismember(mutPositions,ampPositions{o2})); % getting positions that are in amplification region
        [ampGenes, idx1] = unique(mutGenes(loc1)); % getting mutated genes that are in amplification region (including genes overlapping for amplified strain/ not amplified); doing unique so that amplified gene is counted once
        mutGenes = [mutGenes; ampGenes];
        mutClass = [mutClass; repmat("AMP", length(ampGenes), 1)];
        mutType = [mutType; repmat("AMP", length(ampGenes), 1)];
        mutPositions = [mutPositions; mutPositions(loc1(idx1))];
        mutFullName = [mutFullName; ampGenes];
        mutFullInfo = [mutFullInfo; cell(length(ampGenes),3)];
        mutStrain = [mutStrain; repmat(o2,length(ampGenes), 1)];
        if noWriteOver~=2 % not skipping plotting
            for o3 = loc1' % plotting for each position
                plot(a,[firstCircle(mutPositions(o3),1),x(mutPositions(o3))], [firstCircle(mutPositions(o3),2), y(mutPositions(o3))],'k-','LineWidth',2)
            end
        end
    end
end


%% Labeling the gene
% Define the radius and center of the labeling circle
radius = radMult+j*12; % Change the radius as desired
center = [0, 0]; % Change the center coordinates as desired

% Getting labels of genes in the label circle
xlab = center(1) + radius * cos(thetaV);
ylab = center(2) + radius * sin(thetaV);
uniqueGene = unique(mutFullName); % getting unique genes from fullname
locations = zeros(length(uniqueGene),1);
for h = 1:length(uniqueGene)
    posInd = find(strcmp(mutFullName,uniqueGene{h}));
    posGene = mutPositions(posInd(1));
    locations(h) = posGene;
end
[sortLocations, I] = sort(locations);
sortLocations = [sortLocations; max(sortLocations)*2];
sortUG = uniqueGene(I);
sk=0;

geneLabels = {};
for h2 = 1:length(sortUG)
    posInd = find(strcmp(mutFullName,sortUG{h2}));
    sortUG{h2} =strrep(sortUG{h2},'_', '\_'); % replacing underscore so it doesn't subscript
    pg = mutPositions(posInd(1));
    if sortLocations(h2)+10000>sortLocations(h2+1)
        % Concatenate genes that are close together
        textG = sprintf('%s &/or %s', sortUG{h2}, sortUG{h2+1});
        if sk~=1
            geneLabels{h2,1} = textG;
            geneLabels{h2,2} = pg;
            if xlab(pg)>0
                geneLabels{h2,3} = rad2deg(thetaV(pg));
                geneLabels{h2,5} = 'left';
            else
                geneLabels{h2,3} = rad2deg(thetaV(pg))+180;
                geneLabels{h2,5} = 'right';
            end
        end
        sk = 1;
    else
        if sk~=1
            geneLabels{h2,1} = sortUG{h2};
            geneLabels{h2,2} = pg;
            if xlab(pg)>0
                geneLabels{h2,3} = rad2deg(thetaV(pg));
                geneLabels{h2,5} = 'left';
            else
                geneLabels{h2,3} = rad2deg(thetaV(pg))+180;
                geneLabels{h2,5} = 'right';
            end
        end
        sk =0;
    end
    geneLabels{h2,4} = 'k';
end

% making the plot space bigger
if noWriteOver~=2
    xlim(a,[min(xlab*1.5) max(xlab*1.5)])
end
%mutStrain = [mutStrain,zeros(length(mutStrain),1),zeros(length(mutStrain),1)];
%% Compiling strain information
strainInfo = table(mutStrain, mutGenes, mutType, mutClass, mutFullInfo,'VariableNames',["strain", "gene", "type", "category", "info"]);
strainInfo = splitvars(strainInfo, 'info', 'NewVariableNames',["mutation","annotation","description"]);

%% plotting radian lines of genes that are in more than 1 strain
% removing IS element from name
% mutGenesISR = {};
% for i15 = 1:length(mutGenes)
%     tempS = split(mutGenes(i15),'_');
%     mutGenesISR = [mutGenesISR, tempS(~contains(tempS,["IS", "ins"]))'];
% end
catMutGenes = categories(categorical(mutGenes));
allGene = [];
numStrainMutGene = zeros(length(catMutGenes),1);
numStrainLOF = zeros(length(catMutGenes),1);
numStrainMOD = zeros(length(catMutGenes),1);
numStrainAMP = zeros(length(catMutGenes),1);
numStrainJC = zeros(length(catMutGenes),1);
numStrainMC = zeros(length(catMutGenes),1);
numStrainPM= zeros(length(catMutGenes),1);
numStrainCoding = zeros(length(catMutGenes),1);
numStrainIntergenic = zeros(length(catMutGenes),1);
mutationPosition = zeros(length(catMutGenes),1);
nMutDrug = zeros(length(catMutGenes),1);
nMutLOF= zeros(length(catMutGenes),1);
nMutMOD = zeros(length(catMutGenes),1);
nMutAMP = zeros(length(catMutGenes),1);
nMutJC = zeros(length(catMutGenes),1);
nMutMC = zeros(length(catMutGenes),1);
nMutPM = zeros(length(catMutGenes),1);
nMutCoding = zeros(length(catMutGenes),1);
nMutIntergenic = zeros(length(catMutGenes),1);
nMutCodeLOF = zeros(length(catMutGenes),1);
nMutCodeMOD =zeros(length(catMutGenes),1);

h3 = 1;
allGeneLabels = {};
for o1 =1:length(catMutGenes)
    commonGind=contains(mutGenes,catMutGenes{o1}); % getting which strains correspond to the gene
    numStrainMutGene(o1)= nnz(unique(mutStrain(commonGind,1))); % getting # of strains w/ mutations in gene
    % Getting number of strains per category 
    geneStrains = mutStrain(commonGind,1); 
    classes = mutClass(commonGind);
    types= mutType(commonGind);
    numStrainLOF(o1) = sum(nnz(unique(geneStrains(contains(classes,'LOF'),1))));
    numStrainMOD(o1) = sum(nnz(unique(geneStrains(contains(classes,'modification'),1))));
    numStrainAMP(o1) = sum(nnz(unique(geneStrains(contains(classes,'AMP'),1))));
    numStrainJC(o1) = sum(nnz(unique(geneStrains(contains(types,'JC'),1))));
    numStrainMC(o1) = sum(nnz(unique(geneStrains(contains(types,'MC'),1))));
    numStrainPM(o1) = sum(nnz(unique(geneStrains(~contains(types,'JC') & ~contains(types,'MC') & ~contains(types,'AMP'),1)))); % 2026-03-31: Corrected to include subtracting AMP
    numStrainCoding(o1) = sum(nnz(unique(geneStrains(contains(types,'coding'),1)))); % 2026-03-31: added in # of coding/intergenic mutations for predicted mutations that breseq has these calls for 
    numStrainIntergenic(o1) = sum(nnz(unique(geneStrains(contains(types,'intergenic'),1))));

    % Getting number of mutations per category (it's possible that same
    % gene has multiple mutations for the same gene)- by removing unique 
    % Added 2026-04-01
    nMutDrug(o1)= nnz(mutStrain(commonGind,1)); % getting # of mutations for the gene 
    nMutLOF(o1) = sum(nnz(geneStrains(contains(classes,'LOF'),1)));
    nMutMOD(o1) = sum(nnz(geneStrains(contains(classes,'modification'),1)));
    nMutAMP(o1) = sum(nnz(geneStrains(contains(classes,'AMP'),1)));
    nMutJC(o1) = sum(nnz(geneStrains(contains(types,'JC'),1)));
    nMutMC(o1) = sum(nnz(geneStrains(contains(types,'MC'),1)));
    nMutPM(o1) = sum(nnz(geneStrains(~contains(types,'JC') & ~contains(types,'MC') & ~contains(types,'AMP'),1))); 
    nMutCoding(o1) = sum(nnz(geneStrains(contains(types,'coding'),1))); 
    nMutIntergenic(o1) = sum(nnz(geneStrains(contains(types,'intergenic'),1)));
    nMutCodeLOF(o1) = sum(nnz(geneStrains(contains(types,'coding') & contains(classes,'LOF'),1)));
    nMutCodeMOD(o1) = sum(nnz(geneStrains(contains(types,'coding') & contains(classes,'modification'),1)));

    if length(unique(mutStrain(commonGind,1)))>1 % more than 1 strain has mutations in the gene
        allmPos  = mutPositions(commonGind);
        mutationPosition(o1) = allmPos(1); % getting mutation position
        if noWriteOver~=2 % not skipping plotting
            plot(a,[firstCircle(allmPos(1),1),x(allmPos(1))], [firstCircle(allmPos(1),2), y(allmPos(1))],'k-','LineWidth',2)
        end
        % highlighting mutated genes in all strains
        if (length(unique(mutStrain(commonGind,1))))>=length(colInds)
            allPos = mutPositions(contains(mutGenes,catMutGenes{o1}));
            textRot = rad2deg(thetaV(allPos(1)));
            text1 = strrep(catMutGenes{o1}, '_', '\_');% replacing underscore so it doesn't subscript
            allGeneLabels{h3, 1} = text1;
            allGeneLabels{h3, 2} = allmPos(1);
            allGeneLabels{h3, 4} = 'r';
            if xlab(allPos(1))>0
                allGeneLabels{h3, 3} = textRot;
                allGeneLabels{h3,5} = 'left';
            else
                allGeneLabels{h3, 3} = textRot+180;
                allGeneLabels{h3,5} = 'right';
            end
            allGene =[allGene,catMutGenes{o1}];
        end
        h3 = h3+1;
    else
        temp = mutPositions(commonGind);
        mutationPosition(o1) = temp(1);
    end
end

% Plotting gene labels
geneLabels(cellfun(@(x) isempty(x), geneLabels(:,2)),:)=[]; % removing empty rows
if ~isempty(allGeneLabels)
allGeneLabels(cellfun(@(x) isempty(x), allGeneLabels(:,2)),:)=[]; % removing empty rows
try %  removing geneLabels that are in all strains and overlap all strains
    [~,ia] = intersect(cell2mat(geneLabels(:,2)),cell2mat(allGeneLabels(:,2)));
    geneLabels(ia,:) = [];
catch
end
end 
geneLabels(cellfun('isempty', geneLabels(:,1)),:) = [];
for tl = 1:height(geneLabels)
    text(a,xlab(geneLabels{tl,2}),ylab(geneLabels{tl,2}),geneLabels{tl,1},'Rotation',geneLabels{tl,3}, 'HorizontalAlignment',geneLabels{tl,5},'FontSize',14);
end

if ~isempty(allGeneLabels)
    allGeneLabels(cellfun('isempty', allGeneLabels(:,1)),:) = [];
    % if labels are overlapping
    locations = cell2mat(allGeneLabels(:,2));
    labels = allGeneLabels(:,1);
    if isscalar(locations)
        uniqueClusters = locations;
        T = locations; 
    else
        % Perform heirarchical clustering 
        Z = linkage(locations, 'single', 'euclidean');
        % Create clusters based on 10,0000 threshold 
        T = cluster(Z, 'cutoff', 10000,'criterion','distance');
        uniqueClusters = unique(T); 
    end
    for tl = uniqueClusters'
        idx = find(T==tl); 
        if length(labels(T==tl))>1
        textG = strjoin(labels(T==tl), ' and/or ');
        else 
            textG = labels(T==tl);
        end 
        text(a,xlab(allGeneLabels{idx(1),2}),ylab(allGeneLabels{idx(1),2}),textG,'Rotation',allGeneLabels{idx(1),3}, 'HorizontalAlignment',allGeneLabels{idx(1),5},'FontSize',14, 'Color', 'r');
    end 
end


% Compiling all the mutated gene information
totalMutGenes = table(catMutGenes, numStrainMutGene,numStrainLOF, numStrainMOD, numStrainAMP,numStrainJC, numStrainMC, numStrainPM,numStrainCoding, numStrainIntergenic, mutationPosition, ...
    nMutDrug,nMutLOF,nMutMOD,nMutAMP,nMutJC,nMutMC,nMutPM,nMutCoding,nMutIntergenic,nMutCodeLOF,nMutCodeMOD,...
    'VariableNames',{'gene', 'NumStrainsMutated', 'NumStrainsLOF', 'NumStrainsModification','NumStrainsAmplified', 'NumStrainsJC', 'NumStrainsMC', 'NumStrainsPredicted','NumStrainsCoding', 'NumStrainsIntergenic','Position' ...
    'NumMutations', 'MutationsLOF', 'MutationsModification','MutationsAmplified', 'MutationsJC', 'MutationsMC', 'MutationsPM', 'MutationsCoding', 'MutationsIntergenic', 'MutationsCodingLOF', 'MutationsCodingModification'});
totalMutGenes = addvars(totalMutGenes, zeros(height(totalMutGenes),1)+length(colInds), 'NewVariableNames','TotalStrains');

if noWriteOver~=2 % not skipping plotting
    %legend
    strainNames = compareData.Properties.VariableNames(colInds);
    strainNames = strrep(strainNames, '_', ' ');
    % adding in legend
    p=[];
    p(1) = scatter(a,nan, nan,100, 'k','filled','MarkerEdgeColor',[150/255 150/255 150/255]); % LOF -predicted
    p(2) = scatter(a,nan, nan,100, zeros(1,3)+150/255,'filled','MarkerEdgeColor',[150/255 150/255 150/255]); % modification - predicted
    p(3) = scatter(a,nan, nan,100, ones(1,3),'filled','MarkerEdgeColor','k'); % no impact - predicted
    p(4) = plot(a,nan, nan, '*', 'MarkerSize',10, 'LineWidth',2, 'Color', 'k'); % junction
    p(5) = plot(a,nan, nan, '*', 'MarkerSize',10, 'LineWidth',2, 'Color',[150/255 150/255 150/255]); % junction
    p(6)= plot(a,nan, nan, 'square', 'MarkerSize', 10,'LineWidth',2,'Color', 'k'); % missing coverage
    p(7) = plot(a,nan, nan, "_", 'MarkerSize', 10,'LineWidth',2,'Color', 'b'); % amplification
    p(8) = plot(a,nan, nan, "|", 'MarkerSize', 10,'LineWidth',2,'Color', 'b'); % amplification
    legend([legplots, p], [strainNames, "predicted (LOF)", "predicted (modification)", "predicted (no impact)","junction (LOF)", "junction (modification)", "missing coverage", "long amplification", "short amplification"], 'Location', 'southwest')
end

% making COMPARE table for drug: remove unecessary columns
mutIndexA = find(strcmp(compareData.Properties.VariableNames, 'annotation'));
mutIndexG = find(strcmp(compareData.Properties.VariableNames, 'gene'));
compareTableFull = compareData(unique(cell2mat(AstrainInd)),[1:3, colInds, mutIndexA:width(compareData)]);
compareTable = compareData(unique(cell2mat(AstrainInd)),[1:3, colInds, mutIndexG]);

%% saving common mutated genes
drugGenes = {allGene, compareTableFull, drugJC, drugMC, drugAMP, totalMutGenes, strainInfo}; % gene name, # of strains with mutated gene, table(each evolved strain), info for each strain
%%
if noWriteOver~=2
    % Making Amplification table for drug
    drugAMP = drugAMP(~cellfun('isempty',drugAMP)); %removing empty wells
    g2=1;
    temp = cellfun(@(x) x(:,3), drugAMP);
    temp2 = cell2mat(vertcat(temp));
    sNames = cell(height(temp2),1); % intializing
    for g1 = 1:length(drugAMP)
        sNames(g2:(height(drugAMP{g1}{3})+g2-1)) = drugAMP{g1}(1);
        g2 = g2+height(drugAMP{g1}{3});
    end

    if ~isempty(temp2)
        amplificationTable = table(sNames, temp2(:,1), temp2(:,2));
        amplificationTable = renamevars(amplificationTable, ["sNames", "Var2", "Var3"], ["Amplification_Strain", "Start", "End"]);
    else
        amplificationTable =[];
    end

    % Save Tables into Excel file
    xlsFName = sprintf('./BreseqInformation/%s_BreseqResults.xlsx', drugName);
    writetable(compareTable, xlsFName, 'WriteVariableNames',true, 'WriteMode','replacefile');
    if ~isempty(drugJC)
        writetable(drugJC,xlsFName, 'WriteMode','Append','WriteVariableNames',true);
    end
    if ~isempty(drugMC)
        writetable(drugMC,xlsFName, 'WriteMode','Append','WriteVariableNames',true);
    end
    if ~isempty(amplificationTable) % not sure why it's giving a warning- deal with it when it becomes an error
        writetable(amplificationTable,xlsFName, 'WriteMode','Append','WriteVariableNames',true);
    end
    writetable(totalMutGenes, xlsFName, 'WriteMode','append', 'WriteVariableNames',true, 'Sheet', 'TotalMutatedGenes')

    % Plotting Total mutated genes table for drug
    uitable(f,"Data",totalMutGenes,"Position", [Ssize(3)*0.4 0 Ssize(3)*0.5 Ssize(4)],'ColumnWidth', 'auto', 'Units','normalized');

    % saving last figure
    tic
    drawnow; % updating so it will save most up to date
    if noWriteOver ==0 % 0 if you don't want to write over existing file
        if exist(sprintf("./CircaPlots/%s_circa.pdf", drugName), 'file')==0 % Skip if already created
            exportgraphics(a,sprintf("./CircaPlots/%s_circa.pdf", drugName), 'ContentType','vector') % saving circa
            exportapp(f,sprintf("./CircaInformation/%s.pdf", drugName)); %PDF is higher quality
            savefig(f,sprintf("./CircaFigure/%s.fig", drugName), 'compact')
        end
    elseif noWriteOver ==3 % only save the circa Information plots
        exportapp(f,sprintf("./CircaInformation/%s.pdf", drugName)); %PDF is higher quality
    else % if you want to write over exisiting file
        exportgraphics(a,sprintf("./CircaPlots/%s_circa.pdf", drugName), 'ContentType','vector') % saving circa
        exportapp(f,sprintf("./CircaInformation/%s.pdf", drugName)); %PDF is higher quality
        savefig(f,sprintf("./CircaFigure/%s.fig", drugName), 'compact')
    end
    toc
    close all force 
end
fprintf('\n %s analysis done\n',drugName)
end


