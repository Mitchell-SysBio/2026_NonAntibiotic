%% organizingFastQFiles: To move fastq files into folders(SampleName) for sequencing from SeqCenter
% 
%% Initial steps before running
% Add column called "fastQFilename" to sample manifest
% In the column: concatenate
% "OrderNumber_SeqCoastTubeID_"
% You may need to check to make it match the fastq filenames 



%% The code you actually run
% Change the file name and range as needed and only run the second line if you are
% working with a subset of your metadata

sampleMap = readtable("6075_SampleManifest_CarmenLi.csv", 'Range',"E2:R32");
sampleMap = sampleMap(22:25,:); % Run this line if you are only analyzing a subsection of the samples in the metadata

%% Code you don't need to modify
for i = 1:height(sampleMap)
    sampleName = sampleMap.SampleName{i}; % get sample name 
    folderName = sprintf("%s_TubeID%i", sampleName,sampleMap.SeqCoastTubeID(i)); 
    % making folder 
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    fileName = sprintf("%s*",sampleMap.fastQFilename{i});
    movefile(fileName, folderName)
end
