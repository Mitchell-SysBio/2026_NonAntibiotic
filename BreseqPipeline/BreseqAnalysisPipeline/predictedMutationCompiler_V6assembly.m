function [predictedTable] = predictedMutationCompiler_V6assembly(compareFile, gdFile, ancestorName, gdColumnNames)
%% This function gets the predicted mutations from COMPARE file and gets information from the GD file
% Updated 2025-07-09 by Carmen Li
% Correcting error in breseq that makes multiples of 3 LOF 

compareData = readtable(compareFile); % Getting data from COMPARE file 
gdData = readtable(gdFile, "FileType","text",'Delimiter', '\t', 'EmptyValue',Inf); 
% Getting data from TSV file and making empty values Inf because Nan is distinct
% Adding in missing columns in gdData because if statements have been
% hardcoded in by column names
missingColumns = setdiff(gdColumnNames, gdData.Properties.VariableNames);
emptyCol = zeros(height(gdData),length(missingColumns))+Inf;
gdData = addvars(gdData, emptyCol); 
gdData = splitvars(gdData, "emptyCol", "NewVariableNames",missingColumns);
essentialityTable = readtable("Copy-number-and-essentiality-of-all-genes-(E-coli-K-12-substr-MG1655).xlsx");
essentialGenes = essentialityTable{strcmp(essentialityTable.Growth_Status, 'NONE'), 'GeneName'};
%% Combining together the compare Table and gdTable 
% Compare data has additional rows for large deletions because not all the deleted genes fit in the row
% These rows will be removed, but the deleted genes will remain in the
% 'genes_inactivated' column 

%% Group by all columns except the 'drugs' column
% removing variables in order to find unique rows
removeVars(1) = find(strcmp(gdData.Properties.VariableNames, 'title'));
removeVars(2) = find(strcmp(gdData.Properties.VariableNames, 'new_read_count'));
removeVars(3) = find(strcmp(gdData.Properties.VariableNames, 'new_read_count_basis'));
removeVars(4) = find(strcmp(gdData.Properties.VariableNames, 'ref_read_count'));
removeVars(5) = find(strcmp(gdData.Properties.VariableNames, 'ref_read_count_basis'));

gdInds = 1:width(gdData);
gdInds(removeVars) =[];
tempGD = gdData(:,gdInds); % excluding title
[uniqueRows, ~, groupIdx] = unique(tempGD); % Exclude 'drugs' column
drugsGrouped = accumarray(groupIdx, (1:height(tempGD))', [], ...
    @(rows) {gdData.title(rows)});
% add in removing nan rows in compareData to count

if sum(~isnan(compareData.position))~=height(uniqueRows)
    fprintf('ERROR: # of unique rows in GD file is not equal to rows in compareData\n')
else
    fprintf('Correctly matched GD file with compareData\n')
end 
%% Getting gdTable
gdTable = [uniqueRows, table(drugsGrouped, 'VariableNames', {'Samples'})];
gdTable = renamevars(gdTable, "position", "position_gd");

% changing some columns to cells from numeric: because if it's numeric then
% it's empty
classgdTable =  varfun(@class, gdTable, 'OutputFormat', 'cell'); 
colNames = gdTable.Properties.VariableNames;
geneCols = {'genes_inactivated', 'genes_overlapping', 'genes_promoter'}; 
for gc = 1:3
    if strcmp(classgdTable{strcmp(colNames, geneCols{gc})},'double')
        if gc ==1
        gdTable.genes_inactivated = repmat({''},height(gdTable),1);
        elseif gc ==2
        gdTable.genes_overlapping = repmat({''},height(gdTable),1);
        elseif gc ==3
        gdTable.genes_promoter = repmat({''},height(gdTable),1);
        end 
    end 
end 
%% combining gdTable & compareData 
cdInds = nan(height(gdTable),1); % preallocating variable for index of corresponding row in compareData 
classMutation = cell(height(gdTable),1); LOF = cell(height(gdTable),1); MOD = cell(height(gdTable),1); AMP = cell(height(gdTable),1); NoImpact = cell(height(gdTable),1); 
for i = 1:height(gdTable) % need to loop through to assign because can have different mutations at same position
    cdIND = find(compareData{:,"position"}==gdTable{i,"position_gd"});
    mutCat = gdTable{i,"mutation_category"};
    LOFgenes = {}; Modifiedgenes = {}; AMPgenes = {}; noIgenes = {};
    for i1 = cdIND' % looping through each matching index in compareData
        cMutation = compareData{i1,"mutation"}{1}; % getting mutation 
        if contains(mutCat,"snp_")
            % mutation= snp 
            if gdTable{i,"new_seq"}{1} ==cMutation(3) % if its the right new sequence
                cdInds(i) = i1;
                if ~isempty(gdTable{i,"genes_inactivated"}{1})
                   % LOF: if it has inactivated genes
                    classMutation{i} = 'LOF';
                else
                    classMutation{i} = 'modification';
                end 
            end 
        elseif strcmp(mutCat, 'mobile_element_insertion')
            % mutation = IS element by strand IS # (strand +/-) bp change 
            ind = find(cMutation == '(');
            if cMutation(ind+1) == '+'
                temp = 1;
            elseif cMutation(ind+1) == '–'
                temp = -1;
            end 
            if gdTable{i, "strand"} == temp
                cdInds(i) = i1;
                    if ~isempty(gdTable{i,"genes_inactivated"}{1})
                        % LOF: if it has inactivated genes
                        classMutation{i} = 'LOF';
                    elseif rem(gdTable{i,"duplication_size"},3)~=0 && ~isempty(gdTable{i,"genes_promoter"}{1})
                        % LOF has insertion/deletion that isn't mutliple of 3 in promoter region 
                        classMutation{i} = 'LOF';
                    else
                        classMutation{i} = 'modification';
                    end
            end 
        elseif strcmp(mutCat, 'small_indel')
            % mutation = small indel 
            % assigning the right compare table 
            if contains(cMutation, num2str(gdTable{i,"size"}))
                % for if compare table mutation is just a number of bp
                % change
                cdInds(i) = i1;
            elseif contains(cMutation, gdTable{i,"new_seq"}{1}) % if "new_seq" is empty this is also counted in this 
                cdInds(i) = i1;
            else
                fprintf('ERROR: Small Indel is not being accounted for!')
                break 
            end 
                % LOF: if it has inactivated genes or has insertion/
                % deletion that isn't mutliple of 3 in promoter region
                if ~isempty(gdTable{i,"genes_inactivated"}{1})
                    % check if small indel is called correctly
                    if rem(gdTable{i,"size"},3)==0 % Checking if size is multiple of 3 bp
                        % Getting position in coding region
                        temp = strsplit(gdTable{i,"gene_position"}{1}, {'coding (','-'});
                        initialPos = str2num(temp{2});
                        if rem(initialPos,3)==1 % if true then called incorrectly and should be MOD
                            % starting position is 1st position of codon
                             classMutation{i} = 'modification';
                             if isempty(gdTable{i,"genes_overlapping"}{1}) == 1
                                gdTable{i,"genes_overlapping"} = gdTable{i,"genes_inactivated"};
                             else % making comma seperated 
                                 gdTable{i,"genes_overlapping"} = [gdTable{i,"genes_overlapping"}{1}, ',', gdTable{i,"genes_inactivated"}];
                             end 
                             gdTable{i,"genes_inactivated"} = {''};
                        else
                            % LOF: if it has inactivated genes
                            classMutation{i} = 'LOF';
                        end
                    else 
                        % LOF: if it has inactivated genes
                        classMutation{i} = 'LOF';
                    end 
                elseif strcmp(gdTable{i,"type"}{1}, "DEL")
                    if rem(gdTable{i,"size"},3)~=0 && ~isempty(gdTable{i,"genes_promoter"}{1})
                    % LOF: if it has has deletion that isn't mutliple of 3 in promoter region
                    classMutation{i} = 'LOF';
                    else
                        classMutation{i} = 'modification';
                    end
                elseif strcmp(gdTable{i,"type"}{1}, "INS")
                    if ~isempty(gdTable{i,"genes_promoter"}{1}) && rem(gdTable{i,"repeat_new_copies"},3)~=0
                        % LOF: if it has has insertion that isn't mutliple of 3 in promoter region
                        % This one is to account for additional repeat
                        % region
                        classMutation{i} = 'LOF';
                    elseif ~isempty(gdTable{i,"genes_promoter"}{1}) && rem(length(gdTable{i,"new_seq"}{1}),3)~=0
                        % LOF: if it has has insertion that isn't mutliple of 3 in promoter region
                        % This one is for INS that aren't additional
                        % repeats 
                        classMutation{i} = 'LOF';
                    else
                        classMutation{i} = 'modification';
                    end
                else
                    classMutation{i} = 'modification';
                end
        elseif strcmp(mutCat, 'large_deletion') || strcmp(mutCat, 'large_substitution')
            % mutation = large deletion or large substitution
            if contains(cMutation, addComma(gdTable{i,"size"}))
               cdInds(i) = i1;
               if ~isempty(gdTable{i,"genes_inactivated"}{1})
                classMutation{i} = 'LOF';
               else
                    classMutation{i} = 'modification';
               end
            end 
        elseif strcmp(gdTable{i,"type"}{1}, "AMP")
            % mutation = Amplification
            % Will probably needed to be edited when we encounter an AMP
            if contains(cMutation, addComma(gdTable{i,"size"}))
               cdInds(i) = i1;
            end 
            classMutation{i} = 'AMP';
            AMPgenes = [AMPgenes, gdTable{i,"gene_name"}{1}];
        else
            % can't find mutations 
            fprintf('Cannot find mutations in gdTable row: %i\n', i)
        end
        % Classifying genes 
        LOFgenes = gdTable{i,"genes_inactivated"}{1};
        Modifiedgenes = gdTable{i,"genes_overlapping"}{1};

        % Gene classified as promoter only if it did not make
        % genes_inactivated & genes_overlapping list
        if  ~isempty(gdTable{i,"genes_promoter"}{1}) 
            if strcmp(classMutation{i},'LOF')
                LOFgenes = gdTable{i,"genes_promoter"}{1};
            else
                Modifiedgenes = gdTable{i,"genes_promoter"}{1};
            end
        end 

        % Genes classified as no impact if it does not have genes in
        % genes_inactivated, genes_overlapping, genes_promoter
        if isempty(gdTable{i,"genes_inactivated"}{1}) && isempty(gdTable{i,"genes_overlapping"}{1}) && isempty(gdTable{i,"genes_promoter"}{1})
            classMutation{i} = 'noImpact';
            noIgenes = gdTable{i,"gene_name"}{1};
        end 
    end
    % Correcting call for essential gene: making any LOF for essential gene into modification 
    splLOF = strsplit(LOFgenes,',');
    splMOD = strsplit(Modifiedgenes,',');
    % Correcting where name is not truncated for Rhs element protein
    myind = find(contains(splLOF,'rhselementprotein'));
    if ~isempty(myind)
        splLOF(myind) = {splLOF{myind}(end-3:end)};
    end
    eInd = contains(splLOF, essentialGenes); % finding which genes are essential
    if ~isempty(eInd)
        % Recategorizing essential genes to modification 
        LOF(i) = {splLOF(~eInd)};
        if ~isempty(splMOD{1})
            MOD(i) = {[splLOF(eInd), splMOD]};
        else
            MOD(i) = {splLOF(eInd)};
        end 
        if isempty(splLOF(~eInd))
            classMutation{i} = 'modification';
        end 
    else % no essential genes miscatergorized
        LOF(i) = splLOF;
        MOD(i) = splMOD; 
    end 
    AMP(i) = {AMPgenes};
    NoImpact(i) = {noIgenes};
end 

predictedTable = [classMutation compareData(cdInds,:) LOF MOD AMP NoImpact gdTable]; 
colNames = predictedTable.Properties.VariableNames; 
predictedTable = renamevars(predictedTable, colNames(contains(colNames,'Var')), ["MutationClass", "LOFgenes", "Modifiedgenes", "AmplifiedGenes", "noImpactGenes"]);
if predictedTable.position == predictedTable.position_gd
    fprintf("Sucessfully compiled predicted mutations together \n")
else
    fprintf("ERROR mismatching positions\n")
    find((predictedTable.position == predictedTable.position_gd)==0)
end 
%% Removing ancestor mutations, all non NaN
if strcmp(ancestorName, 'no') ~= 1
    ancInds = find(~isnan(predictedTable{:,ancestorName}));
    fprintf("Number of Ancestor mutations removed = %i\n", length(ancInds))
    predictedTable(ancInds, :) = [];
end

% Adding a column at the end: # of strains with mutations in that position 
genes = categorical(predictedTable.gene); % making it into a categorical
totalStrains = cellfun(@(x) length(x), predictedTable{:,"Samples"}); % number of strains per mutation 
predictedTable = addvars(predictedTable, totalStrains); 
% Adding a column at the end: gene names as categorical 
predictedTable = addvars(predictedTable, genes, 'NewVariableNames','catGenes'); 

function numOut = addComma(numIn)
   jf=java.text.DecimalFormat; % comma for thousands, three decimal places
   numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end
end 