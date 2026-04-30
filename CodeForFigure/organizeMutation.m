function [allDrugGeneData, superBugs] = organizeMutation(drugGenesTable,MG1655Position, ABXmech,ABX)
% organizeMutation: takes in the all mutated genes per drug and the position of
% the mutation in MG1655
% Then removes some mutated genes and finds the position of the mutated gene 
% Then organizes it into allDrugGeneData 

% Notes:
% allGenes: Some mutations have no impact if they have no impact then NumStrainsLOF = 0 & NumStrainsModification=0
%% Getting mutated genes + positions in MG1655
drugS = drugGenesTable{:,"DrugName"};
% Removing some mutations that were found in many strains including the controls and is marginal evidence in the ancestor 
noPos = {}; x= {};  
allDrugGeneData = cell(length(drugS),2);
for i = 1:length(drugS)
    % Comparing LFCs 
    drug= drugS(i);
    dgtInd = find(contains(drugGenesTable.DrugName,drug)); % index in drugGenesTable
    mutGenes = drugGenesTable.AllMutatedGeneTable{dgtInd,1}; % getting mutated genes for the drug
    genes = {};genesPositions = []; 
    numStrains = [];skipGenes = []; 
    strainsInterest = [];
    for j =1:height(mutGenes)
        mutatedGene = mutGenes{j,"gene"}{1,1};
        if strcmp(mutatedGene, 'ycfK') 
            skipGenes = [skipGenes, mutatedGene, ', '];
            % skipping because this mutation was found in many strains including the controls and is marginal in the ancestor 
            continue 
        elseif strcmp(mutatedGene, 'tfaE') 
            skipGenes = [skipGenes, mutatedGene, ', '];
            % skipping because this mutation was found in many strains including the controls and is marginal in the ancestor 
            continue 
        elseif strcmp(mutatedGene, 'pinE') 
            skipGenes = [skipGenes, mutatedGene, ', '];
            % skipping because this mutation was found in many strains including the controls and is marginal in the ancestor 
            continue 
        elseif contains(mutatedGene,'gatZ')
            skipGenes = [skipGenes, mutatedGene, ', '];
            % skipping because gatZ is deleted and removal forms a new
            % junction that is captured by breseq
            continue
        elseif contains(mutatedGene,'lacZ')
            skipGenes = [skipGenes, mutatedGene, ', '];
            % skipping because lacZ is deleted and removal forms a new
            % junction that is captured by breseq
            continue
        elseif contains(mutatedGene, 'IS') 
            skipGenes = [skipGenes, mutatedGene, ', '];
            % removing IS genes because it does not have a position
            continue 
        elseif contains(mutatedGene, 'ins')
            skipGenes = [skipGenes, mutatedGene, ', '];
            % removing insertion (ins) genes because it does not have a position
            continue 
        end
        strainsInterest = [strainsInterest, j];
    end 
    numStrains = [mutGenes{strainsInterest, 'NumStrainsMutated'}];
    genes = [mutGenes{strainsInterest,"gene"}]; %removing trailing whitespaces
    mutInterestGenes = mutGenes(strainsInterest,:);
    h = height(mutInterestGenes);
    posGenesS = zeros(h,1);
    posGenesE = zeros(h,1);
    tempStrains = cell(h,1);
    strainInfo = drugGenesTable.StrainInfo{dgtInd,1}; % getting strain info for the drug
    for k = 1:h
            % Getting starting/end positions 
            gene1 = mutInterestGenes{k,"gene"}{1,1};
            if mutInterestGenes{k,"NumStrainsLOF"} == 0 && mutInterestGenes{k,"NumStrainsModification"} == 0
                % Noncoding gene: making the position be the position
                % of the mutation as recorded by breseq
                posGenesS(k) = mutInterestGenes{k,"Position"};
                posGenesE(k) = mutInterestGenes{k,"Position"};
            else
                % some genes have multiple positions because it forms
                % different size proteins so we are just choosing the first
                % one
                posInds = find(strcmp(MG1655Position.Gene,gene1));
                posGenesS(k) = MG1655Position{posInds(1),"Start"};
                posGenesE(k) = MG1655Position{posInds(1),"End"};
            end

            % Getting which strain has the mutation 
            tempStrains{k} = unique(strainInfo.strain(strcmp(strainInfo.gene,gene1)));     
    end

    temp = cell(h,1);
    temp(:,:) = drug;
    ABXmInd = find(contains(ABXmech.Labels,drug));
    % Class = ABX mechanism of the evolved drug; 0 = nonABX inds
    if isempty(ABXmInd)
        class1 = zeros(h,1);
    else
        class1 = zeros(h,1)+ABXmInd;
    end 
    fullTable = addvars(mutInterestGenes, temp, 'Before','gene', 'NewVariableNames','drug');
    fullTable = addvars(fullTable, posGenesS,posGenesE, class1,tempStrains, 'NewVariableNames',{'Start_position', 'End_position', 'class', 'StrainMutated'});
    
    allDrugGeneData{i,1}= fullTable;
    allDrugGeneData(i,2)= drug;
    allDrugGeneData{i,4} = skipGenes; % skipped genes
    if sum(contains(ABX,drug))>0
        allDrugGeneData{i,3}= 1; %ABX
    elseif contains(drug,'Control')
        allDrugGeneData{i,3}= 3; % Controls
    else
        allDrugGeneData{i,3}= 2; %nonABX
    end 
end 

superBugs = {'Chloramphenicol', 'TetracyclineHCl', 'Ciprofloxacin','MitomycinC', 'Cisplatin','PentamidineIsethionate','Dichlorophene','Sertraline'};
end