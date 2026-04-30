function [JC, MC, ancJC,ancMC] = JCMCcompiler_V2(ancestorName, htmlFolder)
%% Set up the Import Options and import the data in order to get the whole csv file
% Updated 2025-01-08 by Carmen Li

% moving into htmlfolder 
curFolder = pwd; 
if ~contains(curFolder,htmlFolder)
    cd(htmlFolder)
end 
dataLines = [1, Inf];

opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "breseqversion0360mutationPredictionsmarginalPredictionssummaryS", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12"];
opts.VariableTypes = ["categorical", "categorical", "string", "string", "string", "string", "string", "double", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName9", "VarName10", "VarName11", "VarName12"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "breseqversion0360mutationPredictionsmarginalPredictionssummaryS", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName9", "VarName10", "VarName11", "VarName12"], "EmptyFieldRule", "auto");

fileNames = dir("*.csv");

% Getting Ancestor data 
csvContents = readtable(sprintf("%s.csv", ancestorName),opts);
% find rows that have new junction evidence
JCind = find(csvContents{:,1}=='Unassigned new junction evidence');
JCevidence = csvContents(JCind+2:end, :);
oddRows = 1:2:height(JCevidence);
evenRows = 2:2:height(JCevidence);
oddSide = JCevidence(oddRows, 3:end);
oddSide.Properties.VariableNames = ["SeqID_A", "position_A", "reads_A","readsCov", "score", "skew","frequency" "annotation_A", "gene_A", "geneDescription_A"];
evenSide = JCevidence(evenRows,2:end-5);
evenSide.Properties.VariableNames = ["SeqID_B", "position_B", "reads_B", "annotation_B", "gene_B", "geneDescription_B"];
ancJC = [oddSide,evenSide];

% find rows that have missing coverage
MCind = find(csvContents{:,1}=='Unassigned missing coverage evidence');
ancMC= csvContents(MCind+2:JCind-1, [4:7 9:end-1]);
ancMC.Properties.VariableNames = ["seq_id", "start", "end", "size", "reads", "gene", "description"];

% Getting JC and Missing Coverage data     
JC = cell(length(fileNames)-1, 4); % intitializing 
MC = cell(length(fileNames)-1, 4); 
%fileNames(strcmp({fileNames.name},sprintf("%s.csv", ancestorName)),:)=[];
for k = 1:length(fileNames)
    file = fileNames(k).name;
    csvContents = readtable(sprintf("%s", file),opts);
    % find rows that have new junction evidence
    JCind = find(csvContents{:,1}=='Unassigned new junction evidence');
    JCevidence = csvContents(JCind+2:end, :);
    oddRows = 1:2:height(JCevidence);
    evenRows = 2:2:height(JCevidence);
    oddSide = JCevidence(oddRows, 3:end);
    oddSide.Properties.VariableNames = ["SeqID_A", "position_A", "reads_A","readsCov", "score", "skew","frequency" "annotation_A", "gene_A", "geneDescription_A"];
    evenSide = JCevidence(evenRows,2:end-5);
    evenSide.Properties.VariableNames = ["SeqID_B", "position_B", "reads_B", "annotation_B", "gene_B", "geneDescription_B"];
    orgJC = [oddSide,evenSide];
    % removing ancestor JC
    for p1 = height(ancJC)
        % ADDED IN: Grabbing SeqID and the position in order to find the
        % right index for both compiled and assembly reference sequences 
        ancInd = find(orgJC.position_A == ancJC.position_A(p1) & orgJC.SeqID_A == ancJC.SeqID_A(p1));
        ancInd2 = find(orgJC.position_B == ancJC.position_B(p1) & orgJC.SeqID_B == ancJC.SeqID_B(p1));
        if ancInd == ancInd2
            orgJC(ancInd,:) = []; % removing ancestor JC
        end
    end 
    %%
    % Determine LOF or modification
    % --inactivating-overlap-fraction = (DEFAULT=0.8)
    % Mutations within this fraction of the length of a gene from its beginning are assigned to the 'genes_inactivating'
    % --promoter-distance <arg>        
    % Mutations upstream and within this distance of the beginning of a gene have it added to their 'genes_promoter' list. (DEFAULT=150)
    annotation = cell(height(orgJC),6); % columns: annotation_A; gene_A; annotation_B; gene_B; LOFgenes; ModifiedGenes
    strMarkers = ['(', '/', ')'];
    orgJC = addvars(orgJC, cell(height(orgJC),1), cell(height(orgJC),1), 'NewVariableNames',{'LOFgenes', 'ModifiedGenes'});
    for p2 = 1:height(orgJC)
        annotation{p2,1} = orgJC{p2,'annotation_A'};
        annotation{p2,2} = orgJC{p2,'annotation_B'};
        annotation{p2,3} = orgJC{p2,'gene_A'};
        annotation{p2,4} = orgJC{p2,'gene_B'};
        LOFgenes ={}; Modifiedgenes = {};
        for p3=1:2
            ant = char(annotation{p2,p3});
            strIndexes =[];
            for p4 = 1:3
                strIndexes = [strIndexes, strfind(ant,strMarkers(p4))];
            end
            
            if contains(ant,'intergenic') 
                % intergenic (+-nt/ +-nt)
                for p5 = 1:3
                    strIndexes = [strIndexes, strfind(ant,strMarkers(p5))];
                end
                for p6 = 1:2
                    if ant(strIndexes(p6)+1) == '‑'
                        geneNames = strsplit(annotation{p2,p3+2}, '/');
                        if str2num(ant(strIndexes(p6)+2:strIndexes(p6+1)-1))<=150
                            LOFgenes = [LOFgenes, geneNames(p6)];
                        else
                            % nothing because the junction is not within
                            % the promoter region so won't affect the gene
                        end
                    end
                end
            else 
                % coding (nt before JC/ total gene size)
                overlapFrac = str2num(ant(strIndexes(1)+1:strIndexes(2)-1))/str2num(ant(strIndexes(2)+1:strIndexes(3)-4));
                if overlapFrac<=0.8
                    LOFgenes = [LOFgenes, annotation{p2,p3+2}];
                else
                    Modifiedgenes = [Modifiedgenes, annotation{p2,p3+2}];
                end
            end
        end
        orgJC{p2,'LOFgenes'} = {LOFgenes};
        orgJC{p2,'ModifiedGenes'} = {Modifiedgenes};
    end 
    
    % Storing JC 
    JC{k,1} = file;
    splitFile = strsplit(file(1:end-4), "_");
    JC{k,2} = splitFile(1);
    JC{k,3} = splitFile(2);
    JC{k,4} = orgJC;

    % find rows that have missing coverage 
    if isempty(JCind) %if no JC
        JCind=height(csvContents) +1;
    end 
    MCind = find(csvContents{:,1}=='Unassigned missing coverage evidence');
    MCevidence = csvContents(MCind+2:JCind-1, [4:7 9:end-1]);
    MCevidence.Properties.VariableNames = ["seq_id", "start", "end", "size", "reads", "gene", "description"];
    % removing ancestor MC by gene because the positions might be slighlty
    % off when barcoding 
    for p = 1:height(ancMC)
        ancInd = find(MCevidence.gene == ancMC.gene(p));
        MCevidence(ancInd,:) = []; % removing ancestor Missing coverage 
    end 

    % Storing MC
    MC{k,1} = file;
    MC{k,2} = splitFile(1);
    MC{k,3} = splitFile(2);
    MC{k,4} = MCevidence;
end 
cd ../ % returning to main folder
end