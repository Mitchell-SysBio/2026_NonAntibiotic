%% gettingBreseqData: To organize the output files from breseq for breseq analysis pipeline 
%% getting index.html from each folder and renaming 
if exist("./Results", "dir") ==0
    mkdir('Results')
    if exist("./Results/htmlfiles", "dir") ==0
        mkdir('Results/htmlfiles')
    end 
end 
folders = dir;
dirFlags = [folders.isdir]; % which is a directory
subFolders = folders(dirFlags); % getting subfolders
subFolderNames = {subFolders(3:end).name}; % skipping . and .. 
for j = 1:length(subFolderNames)
    strainName = subFolderNames{j};
    if strainName~="Results" % checking that it's not the Results folder
        oldFile = sprintf("%s/output/index.html",strainName);
        newFile = sprintf("Results/htmlfiles/%s.html", strainName);
        try
            copyfile(oldFile, newFile) % coping over to Results folder and renaming to strainName
        catch
            continue
        end 
    end
end 
fprintf('Done moving index files \n')
%% getting gd from each folder and renaming it and moving it into results folder for comparing 
if exist("./Results/gdfiles", "dir") ==0
    mkdir('Results/gdfiles')
end 
folders = dir;
dirFlags = [folders.isdir]; % which is a directory
subFolders = folders(dirFlags); % getting subfolders
subFolderNames = {subFolders(3:end).name}; % CHANGE numbering to subfolders you want 
for j = 1:length(subFolderNames)
    strainName = subFolderNames{j};
    if strainName~="Results" % checking that it's not the Results folder
        oldFile = sprintf("%s/output/output.gd",strainName);
        strainType = strsplit(strainName,'_'); 
        newFile = sprintf("Results/gdfiles/%s.gd", strainName);
        try
            copyfile(oldFile, newFile) % coping over to Results folder and renaming to strainName
        catch
            continue 
        end 
    end
end 
fprintf('Done moving gd files \n')
%% getting reference.bam from each folder and renaming 
if exist('./Results/bamfiles', 'dir')==0
    mkdir('Results/bamfiles')
end
folders = dir;
dirFlags = [folders.isdir]; % which is a directory
subFolders = folders(dirFlags); % getting subfolders
subFolderNames = {subFolders(3:end).name}; % CHANGE numbering to subfolders you want 
for j = 1:length(subFolderNames)
    strainName = subFolderNames{j};
    if strainName~="Results" % checking that it's not the Results folder
        oldFile = sprintf("%s/data/reference.bam",strainName);
        newFile = sprintf("Results/bamfiles/%s.bam", strainName);
        try
            copyfile(oldFile, newFile) % coping over to Results folder and renaming to strainName
        catch
            continue
        end 
    end
end 
fprintf('Done moving BAM files \n')