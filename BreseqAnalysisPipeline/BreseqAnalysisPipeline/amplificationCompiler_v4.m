function [amplificationData] = amplificationCompiler_v4(locS, gloS, covCutoff, BAMfolder, ancestorName, ampLengthCutoff, smallCovCutoff)
% Updated 2025-01-08 by Carmen Li
%% Finds location of amplification from BAM file 
% moving into BAMfolder 
curFolder = pwd; 
if ~contains(curFolder,BAMfolder)
    cd(BAMfolder)
end 

files = dir('*.mat');
amplificationData = cell(length(files),3); %intitalizing 
% Getting coverage of ancestor to use for normalization 
iFileAnc = find(contains({files.name}, ancestorName, 'IgnoreCase',true));
data(iFileAnc) = load(files(iFileAnc).name);
wtCov = double(data(iFileAnc).coverage);
wtCov = wtCov/(median(wtCov)); % normalizing by median coverage

for iFile =1:length(files)
    fprintf("%s\n",files(iFile).name)
    data(iFile) = load(files(iFile).name);
    Coverage = double(data(iFile).coverage);
    
    % Getting smoothed coverage
    % adding end to beginning and beginning to end in order to compensate
    % for coverage differences at ends due to circular DNA
    xs1 = Coverage(end-(locS):end); 
    xe2 = Coverage(1:(locS));
    yLocal=smoothdata([xs1 Coverage xe2] ,'movmedian',locS);
    ys1 = Coverage(end-(gloS):end);
    ye2 = Coverage(1:(gloS));
    yGlobal=smoothdata([ys1 Coverage ye2],'movmedian',gloS);
    normCov = log2(yLocal((length(xs1)+1):(end-length(xe2)))./yGlobal((length(ys1)+1:(end-length(ye2))))); % normalize coverage

    ampIndex = find(normCov>covCutoff);% getting index of bases with coverage> 1X median
    
    % Getting large amplification data
    amps = [];
    if ~isempty(ampIndex)
        % saving amplification data 
        fName = files(iFile).name;
        amplificationData{iFile,1} = fName(1:(end-13));
        % Handling if there are multiple large amplifications 
        secAmp = [1 find(diff([ampIndex inf])>1) length(ampIndex)]; % Finding where amplification position is not continuous
        for o5 = 1:length(secAmp)-2
            if o5 == 1
                amps= [amps; ampIndex(secAmp(o5)), ampIndex(secAmp(o5+1))]; % start, end of amplification
            else
                amps= [amps; ampIndex(secAmp(o5)+1), ampIndex(secAmp(o5+1))];
            end
        end
    end 

    % Getting small amplifications that are lost through smoothing 
    nCoverage = Coverage/median(Coverage); % normalizing by median coverage
    wtnCov = nCoverage./wtCov; % dividing by ancestor coverage to remove any spikes in coverage that comes from sequencing method 
    wtnCov(wtnCov == Inf) = 0; % wtCov sometimes has zero
    normCov1 = wtnCov/ median(wtnCov, 'omitnan'); 
    ampInds = find(normCov1>smallCovCutoff);
    
    ampsS = [];
    ampIndsC = []; 
    if ~isempty(ampInds)
        % saving amplification data
        fName = files(iFile).name;
        amplificationData{iFile,1} = fName(1:(end-13));
        % Handling if there are multiple large amplifications
        secAmp1 = [1 find(diff([ampInds inf])>1) length(ampInds)]; % Finding where amplification position is not continuous
        for o5 = 1:length(secAmp1)-2
            if o5 == 1
                ampsS= [ampsS; ampInds(secAmp1(o5)), ampInds(secAmp1(o5+1))]; % start, end of amplification
            else
                ampsS= [ampsS; ampInds(secAmp1(o5)+1), ampInds(secAmp1(o5+1))];
            end
        end
        lenAmps= ampsS(:,2)-ampsS(:,1);
        ampsS(lenAmps<ampLengthCutoff, :)=[]; % removing amplifcations that are shorter than the cutoff
        for i2 = 1:height(ampsS)
            ampIndsC = [ampIndsC, ampsS(i2,1):1:ampsS(i2,2)];
        end 
    end

    allAmp = [amps; ampsS]; % column 1 = start, column 2 = end
    allAmpPositions= [ampIndex, ampIndsC]; % all postions of amplifications 
    fName = files(iFile).name;
    if ~isempty(allAmp)
    % saving amplification data 
        amplificationData{iFile,1} = fName(1:(end-13));
        amplificationData{iFile,2} = allAmpPositions;
        amplificationData{iFile,3} = allAmp;
        
        % Plotting amplification data 
        f=figure('color','w','Position', get(0, 'Screensize'), "Visible", 'off');
        subplot(4,1,1)
        plot(Coverage,'r.')
        xlim([1 length(normCov)])
        hold on
        plot(allAmpPositions, Coverage(allAmpPositions), 'k.')
        title(strrep(fName(1:(end-13)),"_","-"))
        subplot(4,1,2)
        plot(yLocal((length(xs1)+1):(end-length(xe2))),'.')
        hold on
        plot(yGlobal((length(ys1)+1:(end-length(ye2)))),'.')
        xlim([1 length(normCov)])
        xline(allAmp(:,1),'k', 'LineWidth',1.5)
        xline(allAmp(:,2),'k', 'LineWidth',1.5)
        legend('Local', 'Global')
        title(sprintf('Local smoothing = %i; Global smoothing = %i',locS, gloS))
        subplot(4,1,3)
        plot(normCov,'r.')
        xlim([1 length(normCov)])
        hold on
        xline(allAmp(:,1),'k', 'LineWidth',1.5)
        xline(allAmp(:,2),'k', 'LineWidth',1.5)
        yline(covCutoff, 'k', 'LineWidth',1.5)
        title(sprintf('Normalized Coverage cutoff = %d', covCutoff))
        subplot(4,1,4)
        plot(normCov1,'r.')
        xlim([1 length(normCov)])
        hold on
        plot(allAmpPositions, normCov1(allAmpPositions), 'k.')
        yline(smallCovCutoff, 'k', 'LineWidth',1.5)
        title(sprintf('Coverage Normalized by Ancestor and divided by median'))
        saveas(f, sprintf('../AmplificationPlots/%s_Amplification.jpg', fName(1:(end-13))))
    else
        % Saving figure to check even if there is no amplication 
         % Plotting coverage without amplification
        f=figure('color','w','Position', get(0, 'Screensize'), "Visible", 'off');
        subplot(4,1,1)
        plot(Coverage,'r.')
        xlim([1 length(normCov)])
        title(strrep(fName(1:(end-13)),"_","-"))
        subplot(4,1,2)
        plot(yLocal((length(xs1)+1):(end-length(xe2))),'.')
        hold on
        plot(yGlobal((length(ys1)+1:(end-length(ye2)))),'.')
        xlim([1 length(normCov)])
        legend('Local', 'Global')
        title(sprintf('Local smoothing = %i; Global smoothing = %i',locS, gloS))
        subplot(4,1,3)
        plot(normCov,'r.')
        xlim([1 length(normCov)])
        hold on
        yline(covCutoff, 'k', 'LineWidth',1.5)
        title(sprintf('Normalized Coverage cutoff = %d', covCutoff))
        subplot(4,1,4)
        plot(normCov1,'r.')
        xlim([1 length(normCov)])
        hold on
        yline(smallCovCutoff, 'k', 'LineWidth',1.5)
        title(sprintf('Coverage Normalized by Ancestor and divided by median'))
        saveas(f, sprintf('../AmplificationPlots/NoAmplification/%s_Amplification.jpg', fName(1:(end-13))))
        
        % Putting none if there is no amplification 
        amplificationData{iFile,1} = 'none';
    end 
end
cd ../ % returning to main folder
end

