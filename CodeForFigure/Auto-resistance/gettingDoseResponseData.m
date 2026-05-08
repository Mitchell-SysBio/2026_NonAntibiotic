%% gettingDoseResponseData: 
% organizes the dose response data from raw data stored in an excel file
metadataFile = "AllGrowthCurveMetadata_v2.xlsx";

fprintf('Compiling growth curves data \n')
tic
allData = readtable(metadataFile, "Sheet","GrowthCurve");
Data = cell(height(allData),1);
allData = addvars(allData,Data);
allData = addvars(allData,Data, 'NewVariableNames','BlankedData');
allData = addvars(allData,Data, 'NewVariableNames','BlankedMedianData');
drugNames = unique(allData{:, "DrugName"});
strainData  = [];
rows = ["A","B","C","D","E","F","G","H"];
for j = 1:height(allData)
    if ~isempty(allData.SheetName{j})
        j
        sheetData  = readtable("GrowthCurve_All_Results.xlsx", "Sheet", allData.SheetName{j});
        startRow =  allData.StartRow{j};
        endRow =  allData.EndRow{j};
        columns = str2num(allData.Columns{j});
        % getting data using wellnames
        wellNames = [];
        for j1 = find(rows==startRow):find(rows==endRow)
            for j2 = columns
                wellNames = [wellNames,sprintf("%s%i",rows(j1), j2)];
            end
        end

        % removing wells 
        rwell = allData.RemoveWells(j);
        idx = find(strcmp(wellNames, rwell));
        wellNames(idx) = [];

        % Adding Data
        concData = sheetData{:,wellNames};
        allData(j,"Data") = {concData};

        % Adding Normalized Data:
        % 1) Subtract the blank: mean of OD of first 3 readings (0-40mi) for this
        % well in order to compensate for wells that have higher OD intially
        % due to drug cloudiness
        blankedData = concData-mean(concData(1:3,:),'all','omitnan');
        allData(j,"BlankedData") = {blankedData};
        % Adding Median of Replicates of Normalized Data
        allData(j,"BlankedMedianData") = {median(blankedData,2,"omitnan")};

    end
end
toc
allCCDataV1 = allData;
save('Evo_IC50Data_2026-04-22.mat', "allCCDataV1", "allData")
fprintf('Finish compiling growth curves data \n')