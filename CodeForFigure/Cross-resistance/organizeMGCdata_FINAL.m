%% organizeMICdata: 
% to organize maximum growth data from "classifyFullPlates" 
% Organize each evolution experiment into MGC table
% Compile all experiments together
% Get log fold change from WT 

matFiles = ["MGCresultsNonABX_06-Sep-2024.mat", "MGCresultsABX_06-Sep-2024.mat", "MGCresults4evo_05-Sep-2024.mat"];
platemapFiles = ["platemapNonABX.xlsx", "platemapABX.xlsx", "platemap4Evo.xlsx"]; 

%% Organize the data  
MGCtables = cell(3,3); LFCtable = []; maskLFCtable =[];allMGCtable=[];
for i = 1:length(matFiles)
    load(matFiles(i)) % loading mat file
    % Adding in Chlorhexidine data for 4Evo
    if i==3
        load('MGCresultsChlorhexidine_12-Sep-2024.mat')
        %% Adding in new Chlorohexidine results
        maxWTCHL = median(MICtableCHL{4,:}); 
        strainNames = MICtable.Properties.VariableNames;
        chlNames = MICtableCHL.Properties.RowNames;
        CHlvalue = zeros(length(strainNames),1); 
        k2 = 1;
        for k1 = strainNames
            if contains(k1,"WT")
                CHlvalue(k2) = maxWTCHL;
            else
                if k1{1}(2) == '1'
                    varInd = 2;
                elseif k1{1}(2) == '2'
                    varInd = 4;
                elseif k1{1}(2) == '3'
                    varInd = 3;
                elseif k1{1}(2) == '4'
                    varInd = 1;
                end
                rowInd= find(contains(chlNames, [k1{1}(1) '_R' k1{1}(5)]));
                CHlvalue(k2) = MICtableCHL{rowInd,varInd};
            end
            k2= k2+1;
        end 
        MICtable("Chlorhexidine",:) = num2cell(CHlvalue');

        load("MGCresults5FU_RFB_26-Sep-2024.mat")
        %% Adding in new 5-FU results
        maxWTFL = median([MICtable5FU_RFB{["Fluorouracil_A_R1", "Fluorouracil_B_R1", "Fluorouracil_C_R1", "Fluorouracil_D_R1", "Fluorouracil_A_R2", "Fluorouracil_B_R2"],"ANC.jpg"}; ...
            MICtable5FU_RFB{["Fluorouracil_A_R1", "Fluorouracil_B_R1", "Fluorouracil_C_R1", "Fluorouracil_D_R1"], "ANC_2.jpg"}]); 
        strainNames = MICtable.Properties.VariableNames;
        FLNames = MICtable5FU_RFB.Properties.RowNames;
        FLvalue = zeros(length(strainNames),1); 
        k2 = 1;
        for k1 = strainNames
            if contains(k1,"WT")
                FLvalue(k2) = maxWTFL;
            else
                if k1{1}(2) == '2'
                    varInd = 6; % Streptozotocin
                    rowInd= find(contains(FLNames, ['Fluorouracil_' k1{1}(1) '_R' k1{1}(5)]));
                    FLvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                elseif k1{1}(2) == '4'
                    varInd = 3; % Ciprofloxacin
                    rowInd= find(contains(FLNames, ['Fluorouracil_' k1{1}(1) '_R' k1{1}(5)]));
                    FLvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                elseif k1{1}(2) == '1'
                    varInd = 4; % Lamotrigine
                    rowInd= find(contains(FLNames, ['Fluorouracil_' k1{1}(1) '_R' k1{1}(5)]));
                    FLvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                elseif k1{1}(2) == '3'
                    varInd = 5; % Sertraline
                    rowInd= find(contains(FLNames, ['Fluorouracil_' k1{1}(1) '_R' k1{1}(5)]));
                    FLvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                else 
                    FLvalue(k2) = nan;
                end
                
            end
            k2= k2+1;
        end 
        MICtable("Fluorouracil",:) = num2cell(FLvalue');
        %% Adding in new Rifabutin results
        maxWTRFB = median([MICtable5FU_RFB{["Rifubutin_A_R1", "Rifubutin_B_R1", "Rifubutin_C_R1", "Rifubutin_D_R1", "Rifubutin_A_R2", "Rifubutin_B_R2"],"ANC.jpg"};...
            MICtable5FU_RFB{["Rifubutin_A_R1", "Rifubutin_B_R1", "Rifubutin_C_R1", "Rifubutin_D_R1"],"ANC_2.jpg"}]); 
        strainNames = MICtable.Properties.VariableNames;
        RFBNames = MICtable5FU_RFB.Properties.RowNames;
        RFBvalue = zeros(length(strainNames),1); 
        k2 = 1;
        for k1 = strainNames
            if contains(k1,"WT")
                RFBvalue(k2) = maxWTRFB;
            else
                if k1{1}(2) == '5'
                    varInd = 2; % Streptozotocin
                    rowInd= find(contains(RFBNames, ['Rifubutin_' k1{1}(1) '_R' k1{1}(5)]));
                    RFBvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                elseif k1{1}(2) == '4'
                    varInd = 1; % Ciprofloxacin
                    rowInd= find(contains(RFBNames, ['Rifubutin_' k1{1}(1) '_R' k1{1}(5)]));
                    RFBvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                elseif k1{1}(2) == '1'
                    varInd = 4; % Lamotrigine
                    rowInd= find(contains(RFBNames, ['Rifubutin_' k1{1}(1) '_R' k1{1}(5)]));
                    RFBvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                 elseif k1{1}(2) == '3'
                    varInd = 5; % Sertraline
                    rowInd= find(contains(RFBNames, ['Rifubutin_' k1{1}(1) '_R' k1{1}(5)]));
                    RFBvalue(k2) = MICtable5FU_RFB{rowInd,varInd};
                else 
                    RFBvalue(k2) = nan;
                end
                
            end
            k2= k2+1;
        end 
        MICtable("Rifubutin",:) = num2cell(RFBvalue');
        %%
    end

    platemapFileName = platemapFiles(i);
    concMap = readmatrix(platemapFileName,'Range','B14:M21');
    drugMap = readcell(platemapFileName,'Range','B3:M10');
    drugMechanism = readcell(platemapFileName,'Range','B36:M43');
    
    fileNames = MICtable.Properties.VariableNames;

    MGCdrugTable = rows2vars(MICtable); % transpose so columns are drugs, rows are plates
    %CalldrugTable = rows2vars(CallTable);
    fileNames = MICtable.Properties.VariableNames;

    repI = [];
    strainName = [];

    % removing ".jpg" or ".png" ending from fileNames
    fileNames = strrep(fileNames,".png","");
    fileNames = strrep(fileNames,".jpg","");

    % getting replicate identifier & strain name of plate
    for o = 1:length(fileNames)
        fileN = fileNames{o};
        if contains(fileN,"WT")
            repI = [repI, str2num(fileN(end))];
            strainName = [strainName, "WT"];
        else
            a = strsplit(fileN,"_");
            %repI = [repI, str2num(a{2})];
            repI = [repI, 1]; % only 1 replicate
            strainName = [strainName, string(a{1})];
        end
    end
    MGCdrugTable.("ReplicateIdentifier") = repI';
    MGCdrugTable.("Strain") = strainName';
    MGCdrugTable = renamevars(MGCdrugTable, "OriginalVariableNames","PlateName");
    plateNames = MGCdrugTable.PlateName;
    plateNames= cellfun(@(x) strrep(x,"_","-"), plateNames);
    MGCdrugTable.PlateName = plateNames;
    MGCdrugTable = movevars(MGCdrugTable,"Strain","After","PlateName");
    MGCdrugTable = movevars(MGCdrugTable,"ReplicateIdentifier","After","Strain");
    % Importing in evo strain information
    evoMap = readcell(platemapFileName,'Range','B25:M32');
    evoMechanisms = readcell(platemapFileName,'Range','B47:M54');
    evoDrugs = cell(length(strainName),1); 
    evoMechs = cell(length(strainName),1); 
    evoType = zeros(length(strainName),1); %1 = ABX, 2 = non-ABX, 3 = control, 0 = WT
    Alphabet = 'ABCDEFGH';
    for sN = 1:length(strainName)
        xsN= strainName{sN};
        evoDrug = evoMap(find(xsN(1)==Alphabet), str2num(xsN(2:end)));
        evoMech = evoMechanisms(find(xsN(1)==Alphabet), str2num(xsN(2:end)));
        if isempty(evoDrug)
            evoDrugs{sN,1} = "WT";
            evoMechs{sN,1} = "WT"; 
        else
            evoDrugs{sN,1} = evoDrug;
            evoMechs{sN,1} = evoMech;
            if i == 1
                evoType(sN,1) = 2; 
            elseif i == 2
                evoType(sN,1) = 1; 
            elseif i ==3
                if contains(evoDrug, "Ciprofloxacin")
                    evoType(sN,1) = 1; 
                else
                    evoType(sN,1) = 2; 
                end 
            end 
            if contains(evoDrug, "Control")
                evoType(sN,1) = 3;
            end 
        end
    end
    MGCdrugTable.("EvolvedDrug")= string(evoDrugs);
    MGCdrugTable.("EvolvedDrugMechanism") = evoMechs;
    MGCdrugTable.("Type")= evoType; %1 = ABX, 2 = non-ABX, 3 = control, 0 = WT
    MGCdrugTable.("BMDSet") = zeros(length(evoType),1)+i; %1 = nonABX, 2 = ABX, 3 = 4Evo
    MGCdrugTable = movevars(MGCdrugTable,["EvolvedDrug", "EvolvedDrugMechanism", "Type", "BMDSet"],"After","ReplicateIdentifier");
    drugStart = 8; % column at which ABX MGC starts 
    %% Get log fold change 
    % choosing the mode aka lowest MGC for the replicates 
    % getting the fold of dilutions per drug
    
    % getting concentration of BMD plate: using information from
    % "platemapNonABX.xlsx" because it's the same for all BMD plates 
    platemapFileName = platemapFiles(1);
    concMap = readmatrix(platemapFileName,'Range','B14:M21');
    drugMap = readcell(platemapFileName,'Range','B3:M10');
    sortConcMap = [sort(concMap(1:4,:)), sort(concMap(5:8,:))]';
    foldDil = table;
    foldDil.Drug = [drugMap(1,:), drugMap(5,:)]';
    foldDil.DilutionFactor = (sortConcMap(:,2)./sortConcMap(:,1)); % dilution factor for each drug
    foldDil.Threshold = max(sortConcMap,[],2)./min(sortConcMap,[],2);
    foldDil.Concentrations = sortConcMap;
    
    MGCdrugAlt = removevars(MGCdrugTable, "PositiveControl");
    drugList = MGCdrugAlt.Properties.VariableNames; 
    drugList(1:7) = [];
    wtInd = find(MGCdrugTable.Strain == "WT");
    wtData = MGCdrugAlt{wtInd,drugStart:end};
    rifInd = find(contains(drugList, "Rifampicin"));
    
    % using median because there is a lot of variability in the WT
    % Though I am using mode for N=3 per sample but median and mode should be
    % pretty much the same then 
    rifMed = median(wtData(:,rifInd), 'all'); % finding median of both replicates of rifampicin for all
    medianData = median(wtData, 1, 'omitnan'); % returning median of each column
    medianData(rifInd(1)) = rifMed;
    medianData(rifInd(2)) = [];
    maxWTdata = medianData; % using the median(maximum growth) of all replicates for wild type to get concensus across replicates  
    
    MGCdrugAlt(wtInd, :) = []; %removing wt 
    
    drugOrder = [drugMap(1,:), drugMap(5,:)];
    drugMechOr = [drugMechanism(1,:), drugMechanism(5,:)];
    removeInds = [find(drugOrder == "PositiveControl")];
    uniStrain = unique(string(MGCdrugAlt.EvolvedDrug));
    try
    foldDil(removeInds,:)=[];
    catch
        % nothing because can't removeInds if it's already removed 
    end 
    foldDil = sortrows(foldDil,"Drug","ascend");
    maxGdata = MGCdrugAlt{:,drugStart:end};
    cleanGdata = zeros(length(uniStrain),width(maxGdata)-1);
    evoDrugChar = cell(length(uniStrain),1); 
    evoDrugMech = cell(length(uniStrain),1); 
    evoDrug = cell(length(uniStrain),1); 
    evoType = zeros(length(uniStrain),1);
    for k = 1:length(uniStrain)
        evoDrugChar{k}= uniStrain{k};
        evoDrug{k} = uniStrain{k}(1:end-3);
        ind = find(contains(MGCdrugAlt.EvolvedDrug,uniStrain(k)));
        evoDrugMech{k} = MGCdrugAlt{ind(1), "EvolvedDrugMechanism"}{1};
        evoType(k) = MGCdrugAlt{ind(1), "Type"};
        strainData = maxGdata(ind,:);
        
        % returns "majority" vote, if all 3 are different then the lowest concentrations
        rifInd = find(contains(drugList, "Rifampicin"));
        rifMode = mode(strainData(:,rifInd), 'all'); % finding mode of both replicates of rifampicin for all 
        modeData = mode(strainData, 1); % returning mode of each column 
        modeData(rifInd(1)) = rifMode; 
        modeData(rifInd(2)) = []; 
    
        cleanGdata(k,:) = modeData;
   end 
    
    datDiv = (cleanGdata./maxWTdata); % divide by maximum wt growth  
    logData = log2(datDiv);
    drugNames = MGCdrugAlt.Properties.VariableNames([drugStart:25 27:end]);
    drugList(rifInd(2)) = [];
    [sortDL sortInd] = sort(drugList);
    LFC = array2table(logData(:,sortInd), 'VariableNames', sortDL, 'RowNames',evoDrugChar');
    LFC.Drug = evoDrug;
    LFC.Mechanism = evoDrugMech; 
    LFC.Type = evoType; 
    LFC = movevars(LFC, ["Type", "Fluorouracil"], "After","Mechanism");
    
    % making LFC table with masking marginal LFC< 2 fold 
    maskLogData = logData;
    maskLogData(maskLogData<1 & maskLogData>-1) = 0;
    maskLFC = array2table(maskLogData(:,sortInd), 'VariableNames', sortDL, 'RowNames',evoDrugChar');
    maskLFC.Drug = evoDrug;
    maskLFC.Mechanism = evoDrugMech; 
    maskLFC.Type = evoType; 
    maskLFC = movevars(maskLFC, ["Fluorouracil"], "After","Mechanism");

    MGCtables{i,1} = matFiles(i); % mat file name
    MGCtables{i,2} = MGCdrugTable; % organized MGC table
    MGCtables{i,3} = MICtable; % orginal MGC table from "classifyFullPlates.m"
    MGCtables{i,4} = LFC; % Log2 fold change from WT
    MGCtables{i,5} = maskLFC; % Masked log2 fold change from WT
    MGCtables{i,6} = {maxWTdata(sortInd),sortDL}; % MGC of WT used for LFC 
    % Compiling together all MGCtables
    allMGCtable = [allMGCtable; MGCdrugTable];
    LFCtable = [LFCtable; LFC];
    maskLFCtable = [maskLFCtable; maskLFC];
end

%% Getting mechanism colors 
% ABXmech: Labels = BMD drug names; Colors = correspond to mechanism; Mechanism = ABX mechanism 
mechColor = [186,85,211; 125,25,207; 220,20,60; 255,140,0; 255,215,0; 0,100,100; 0,128,0; 0,0,180; 240,128,128; 218,112,214; 100,100,100; 0,0,0; 234,155,74; 73,88,157] ;
mechColor = mechColor./255;
mechCname = ["30s","50s","DNA","broad","cellwall","fattyacid","folicacid", "membrane", "nitrofuran", "protein", "unknown", "chemotherapy", "RNA", "cellmembrane"]; 
drugNames = maskLFCtable.Properties.VariableNames([1:21]);
drugMechColor = {}; mechs = {}; 
for l = 1:length(drugNames)
    myind = find(contains(drugOrder,drugNames(l)));
    drugMechColor{l} = mechColor(drugMechOr{myind(1)}==mechCname,:);
    mechs{l} = drugMechOr{myind(1)};
end 
ABXmech.Labels = drugNames;
ABXmech.Colors = drugMechColor;
ABXmech.Mechanism = mechs;


%% Removing Carmofur from "allLFC_092624.mat"
% Not removing from MGCtables 
height(allMGCtable)
height(LFCtable)
height(maskLFCtable)
carNames = ["Carmofur  L1", "Carmofur  L2", "Carmofur  L3", "Carmofur  L4"];
for c = 1:length(carNames)
    ind = strcmp(allMGCtable.EvolvedDrug, carNames(c));
    allMGCtable(ind,:) = [];
    ind2 = strcmp(LFCtable.Properties.RowNames, carNames(c));
    LFCtable(ind2,:) = [];
    ind2 = strcmp(maskLFCtable.Properties.RowNames, carNames(c));
    maskLFCtable(ind2,:) =[];
end 
height(allMGCtable)
height(LFCtable)
height(maskLFCtable)
%% Organize the data 

matFiles = ["MGCresultsCarmofur_02-Jul-2025"];
platemapFiles = ["platemapCarmofurRepeat.xlsx"]; 

for i = 1
    load(matFiles(i)) % loading mat file
    platemapFileName = platemapFiles(i);
    concMap = readmatrix(platemapFileName,'Range','B14:M21');
    drugMap = readcell(platemapFileName,'Range','B3:M10');
    drugMechanism = readcell(platemapFileName,'Range','B36:M43');
    
    fileNames = MICtable.Properties.VariableNames;

    MGCdrugTable = rows2vars(MICtable); % transpose so columns are drugs, rows are plates
    %CalldrugTable = rows2vars(CallTable);
    fileNames = MICtable.Properties.VariableNames;

    repI = [];
    strainName = [];

    % removing ".jpg" or ".png" ending from fileNames
    fileNames = strrep(fileNames,".png","");
    fileNames = strrep(fileNames,".jpg","");

    % getting replicate identifier & strain name of plate
    for o = 1:length(fileNames)
        fileN = fileNames{o};
        if contains(fileN,"WT")
            repI = [repI, str2num(fileN(end))];
            strainName = [strainName, "WT"];
        else
            a = strsplit(fileN,"_");
            %repI = [repI, str2num(a{2})];
            repI = [repI, 1]; % only 1 replicate
            strainName = [strainName, string(a{1})];
        end
    end
    MGCdrugTable.("ReplicateIdentifier") = repI';
    MGCdrugTable.("Strain") = strainName';
    MGCdrugTable = renamevars(MGCdrugTable, "OriginalVariableNames","PlateName");
    plateNames = MGCdrugTable.PlateName;
    plateNames= cellfun(@(x) strrep(x,"_","-"), plateNames);
    MGCdrugTable.PlateName = plateNames;
    MGCdrugTable = movevars(MGCdrugTable,"Strain","After","PlateName");
    MGCdrugTable = movevars(MGCdrugTable,"ReplicateIdentifier","After","Strain");
    % Importing in evo strain information
    evoMap = readcell(platemapFileName,'Range','B25:M32');
    evoMechanisms = readcell(platemapFileName,'Range','B47:M54');
    evoDrugs = cell(length(strainName),1); 
    evoMechs = cell(length(strainName),1); 
    evoType = zeros(length(strainName),1); %1 = ABX, 2 = non-ABX, 3 = control, 0 = WT
    Alphabet = 'ABCDEFGH';
    for sN = 1:length(strainName)
        xsN= strainName{sN};
        evoDrug = evoMap(find(xsN(1)==Alphabet), str2num(xsN(2:end)));
        evoMech = evoMechanisms(find(xsN(1)==Alphabet), str2num(xsN(2:end)));
        if isempty(evoDrug)
            evoDrugs{sN,1} = "WT";
            evoMechs{sN,1} = "WT"; 
        else
            evoDrugs{sN,1} = evoDrug;
            evoMechs{sN,1} = evoMech;
            if i == 1
                evoType(sN,1) = 2; 
            elseif i == 2
                evoType(sN,1) = 1; 
            elseif i ==3
                if contains(evoDrug, "Ciprofloxacin")
                    evoType(sN,1) = 1; 
                else
                    evoType(sN,1) = 2; 
                end 
            end 
            if contains(evoDrug, "Control")
                evoType(sN,1) = 3;
            end 
        end
    end
    MGCdrugTable.("EvolvedDrug")= string(evoDrugs);
    MGCdrugTable.("EvolvedDrugMechanism") = evoMechs;
    MGCdrugTable.("Type")= evoType; %1 = ABX, 2 = non-ABX, 3 = control, 0 = WT
    MGCdrugTable.("BMDSet") = zeros(length(evoType),1)+i; %1 = nonABX, 2 = ABX, 3 = 4Evo
    MGCdrugTable = movevars(MGCdrugTable,["EvolvedDrug", "EvolvedDrugMechanism", "Type", "BMDSet"],"After","ReplicateIdentifier");
    drugStart = 8; % column at which ABX MGC starts 
    %% Get log fold change 
    % choosing the mode aka lowest MGC for the replicates 
    % getting the fold of dilutions per drug
    
    % getting concentration of BMD plate: using information from
    % "platemapNonABX.xlsx" because it's the same for all BMD plates 
    platemapFileName = platemapFiles(1);
    concMap = readmatrix(platemapFileName,'Range','B14:M21');
    drugMap = readcell(platemapFileName,'Range','B3:M10');
    sortConcMap = [sort(concMap(1:4,:)), sort(concMap(5:8,:))]';
    foldDil = table;
    foldDil.Drug = [drugMap(1,:), drugMap(5,:)]';
    foldDil.DilutionFactor = (sortConcMap(:,2)./sortConcMap(:,1)); % dilution factor for each drug
    foldDil.Threshold = max(sortConcMap,[],2)./min(sortConcMap,[],2);
    foldDil.Concentrations = sortConcMap;
    
    MGCdrugAlt = removevars(MGCdrugTable, "PositiveControl");
    drugList = MGCdrugAlt.Properties.VariableNames; 
    drugList(1:7) = [];
    wtInd = find(MGCdrugTable.Strain == "WT");
    wtData = MGCdrugAlt{wtInd,drugStart:end};
    rifInd = find(contains(drugList, "Rifampicin"));
    
    % using median because there is a lot of variability in the WT
    % Though I am using mode for N=3 per sample but median and mode should be
    % pretty much the same then 
    rifMed = median(wtData(:,rifInd), 'all'); % finding median of both replicates of rifampicin for all
    medianData = median(wtData, 1, 'omitnan'); % returning median of each column
    medianData(rifInd(1)) = rifMed;
    medianData(rifInd(2)) = [];
    maxWTdata = medianData; % using the mode(maximum growth) of all replicates for wild type to get concensus across replicates  
    
    MGCdrugAlt(wtInd, :) = []; %removing wt 
    
    drugOrder = [drugMap(1,:), drugMap(5,:)];
    drugMechOr = [drugMechanism(1,:), drugMechanism(5,:)];
    removeInds = [find(drugOrder == "PositiveControl")];
    uniStrain = unique(string(MGCdrugAlt.EvolvedDrug));
    try
    foldDil(removeInds,:)=[];
    catch
        % nothing because can't removeInds if it's already removed 
    end 
    foldDil = sortrows(foldDil,"Drug","ascend");
    maxGdata = MGCdrugAlt{:,drugStart:end};
    cleanGdata = zeros(length(uniStrain),width(maxGdata)-1);
    evoDrugChar = cell(length(uniStrain),1); 
    evoDrugMech = cell(length(uniStrain),1); 
    evoDrug = cell(length(uniStrain),1); 
    evoType = zeros(length(uniStrain),1);
    for k = 1:length(uniStrain)
        evoDrugChar{k}= uniStrain{k};
        evoDrug{k} = uniStrain{k}(1:end-3);
        ind = find(contains(MGCdrugAlt.EvolvedDrug,uniStrain(k)));
        evoDrugMech{k} = MGCdrugAlt{ind(1), "EvolvedDrugMechanism"}{1};
        evoType(k) = MGCdrugAlt{ind(1), "Type"};
        strainData = maxGdata(ind,:);
        
        % returns "majority" vote, if all 3 are different then the lowest concentrations
        rifInd = find(contains(drugList, "Rifampicin"));
        rifMode = mode(strainData(:,rifInd), 'all'); % finding mode of both replicates of rifampicin for all 
        modeData = mode(strainData, 1); % returning mode of each column 
        modeData(rifInd(1)) = rifMode; 
        modeData(rifInd(2)) = []; 
    
        cleanGdata(k,:) = modeData;
   end 
    
    datDiv = (cleanGdata./maxWTdata); % divide by maximum wt growth  
    logData = log2(datDiv);
    drugNames = MGCdrugAlt.Properties.VariableNames([drugStart:25 27:end]);
    drugList(rifInd(2)) = [];
    [sortDL sortInd] = sort(drugList);
    LFC = array2table(logData(:,sortInd), 'VariableNames', sortDL, 'RowNames',evoDrugChar');
    LFC.Drug = evoDrug;
    LFC.Mechanism = evoDrugMech; 
    LFC.Type = evoType; 
    LFC = movevars(LFC, ["Type", "Fluorouracil"], "After","Mechanism");
    
    % making LFC table with masking marginal LFC< 2 fold 
    maskLogData = logData;
    maskLogData(maskLogData<1 & maskLogData>-1) = 0;
    maskLFC = array2table(maskLogData(:,sortInd), 'VariableNames', sortDL, 'RowNames',evoDrugChar');
    maskLFC.Drug = evoDrug;
    maskLFC.Mechanism = evoDrugMech; 
    maskLFC.Type = evoType; 
    maskLFC = movevars(maskLFC, ["Fluorouracil"], "After","Mechanism");

    MGCtables{i,1} = matFiles(i); % mat file name
    MGCtables{i,2} = MGCdrugTable; % organized MGC table
    MGCtables{i,3} = MICtable; % orginal MGC table from "classifyFullPlates.m"
    MGCtables{i,4} = LFC; % Log2 fold change from WT
    MGCtables{i,5} = maskLFC; % Masked log2 fold change from WT
    MGCtables{i,6} = {maxWTdata(sortInd),sortDL}; % MGC of WT used for LFC 
    % Compiling together all MGCtables
    allMGCtable = [allMGCtable; MGCdrugTable];
    LFCtable = [LFCtable; LFC];
    maskLFCtable = [maskLFCtable; maskLFC];
end

%% saving variables
allMGCtable = renamevars(allMGCtable, 'Rifubutin', 'Rifabutin');
LFCtable = renamevars(LFCtable, 'Rifubutin', 'Rifabutin');
maskLFCtable = renamevars(maskLFCtable, 'Rifubutin', 'Rifabutin');
save("allLFC_021926.mat", "allMGCtable", "LFCtable", "maskLFCtable", "ABXmech", "MGCtables")
fprintf('saved!\n')
