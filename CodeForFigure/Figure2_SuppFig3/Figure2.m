%% Fig2 & Supplementary Figure 3-4: to graph heteroresistance
% Fig 2A: bar and radar plots of rifampicin
% Fig 2B: resistance profile radar plots of all drugs 
% Supplementary Figure 3: radar plots of controls 
% Supplementary Figure 4: sensitivity profile radar plots of all drugs 
%% Load data
cd ../
load('allLFC_021926.mat')
load('ABXmechColorMap.mat')
cd Figure2_SuppFig3-4/ 
% "allMGCtable" = col 1: mat file name, col 2 : organized MGC table, col 3:orginal MGC table from "classifyFullPlates.m"
% "LFCtable" = table of log2 fold change in comparision to the wt 
% "maskLFCtable" = table of log2 fold change in comparision to the wt, with
% masking marginal LFC< 2 fold 
% ABXmech: Labels = BMD drug names; Colors = correspond to mechanism; Mechanism = ABX mechanism 

uniMech = unique(ABXmech.Mechanism);
rearrangeI = [2,1,8,3,7,4,5,6];
uniMech = uniMech(rearrangeI); % reorganing so it's in desired order 
%% Fig 2A: Radar plots for abx + bar for each replicate
% get # of cross-resistant ABX for each category 
uniDrugEvo = {"Rifampicin"}; 
mechCrossRes = zeros(length(uniDrugEvo),8); % getting # of cross resistance per mechanism 
mechColorAll = cell(length(uniMech),1); mechInds=cell(length(uniMech),1);
for k = 1:length(uniMech)
    mechInds{k} = ABXmech.Labels((contains(ABXmech.Mechanism, uniMech{k})));
    mechColorAll{k} = ABXmech.Colors((contains(ABXmech.Mechanism, uniMech{k})));
     for g = 1:length(uniDrugEvo)
             mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),mechInds{k}};
             mechCrossRes(g,k) = sum((mydata>0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
     end
end

mechNames= {sprintf("50s\n(3 drugs)"),sprintf("30s\n(2 drugs)"),sprintf("nitrofuran\n(2 drugs)"),sprintf("DNA\n(3 drugs)"),sprintf("folic acid\n(3 drugs)"),sprintf("RNA\n(2 drugs)"),sprintf('cell membrane\n (3 drugs)'),sprintf("cell wall\n(3 drugs)")}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);


fig = figure('color','w'); 
tiledlayout(4,1)
reorderABX = horzcat(mechInds{:});
reorderColors = horzcat(mechColorAll{:});
for s1=1:4
    nexttile
    hold on;
    mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),reorderABX};
    for s2 = 1:length(reorderABX)
        bar(reorderABX(s2),mydata(s1,s2), 'FaceColor',reorderColors{s2})
    end
    ylim([-3 3])
    grid on;box on;
    title(sprintf('Strain %i', s1))
    ylabel('MGC increase (log2 fold change)')
end 
exportgraphics(fig,'Fig2A.svg','ContentType', 'vector');

fig = figure;
for g = 1:length(uniDrugEvo)
    data = mechCrossRes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    
    % spider plot
    spider_plot(data, 'AxesLabels', mechNames, 'AxesDisplay', 'one','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 12,'LabelFont', 'Avenir','LabelFontSize', 10,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,uniDrugEvo{g}));

        title(sprintf("%s", uniDrugEvo{g}),'FontSize',16, 'FontName','Avenir','Color',ABXmech.Colors{tInd})
end
exportgraphics(fig,'Fig2B.svg','ContentType', 'vector');
%% Figure 2B: ABX evo Strains: Radar plots per evoDrug: categorized by ABX mechanism 
% get # of cross-resistant ABX for each category 
% if want to start 0 in middle change axes offset to 0
abxLFC = maskLFCtable(maskLFCtable.Type==1,:); 

[uniDrugEvo, uniInds] = unique(abxLFC.Drug(abxLFC.Type==1));
mechUni = string(abxLFC.Mechanism(uniInds));
desireMechOrder = ["50s", "30s", "nitrofuran", "DNA", "folicacid","RNA","cellmembrane", "cellwall"];
orgInds = [];
for j9=1:length(desireMechOrder) 
    orgInds = [orgInds; find(contains(mechUni, desireMechOrder(j9)))];
end
uniDrugEvo = uniDrugEvo(orgInds); % REORDER via mechanism in same order as radar plot
mechCrossRes = zeros(length(uniDrugEvo),8); % getting # of cross resistance permechanism 
for k = 1:length(uniMech)
    for g = 1:length(uniDrugEvo)
        mydata = abxLFC{(contains(abxLFC.Drug,uniDrugEvo{g})),mechInds{k}};
        mechCrossRes(g,k) = sum((mydata>0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
    end
end

mechNames= {"50s","30s","nitrofuran","DNA","folic acid","RNA",sprintf('cell      \n membrane'),"cell wall"}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);
fig = figure('color','w','Position', get(0, 'Screensize')); 
        spI = 1; 
        r=4; NumEvoLines=6;
for g = 1:length(uniDrugEvo)
    data = mechCrossRes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    subplot(4,5,spI)
    
    % spider plot
     spider_plot(data,  'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 6,'LabelFont', 'Avenir Black','LabelFontSize', 8,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,erase(uniDrugEvo{g}, " ")));
    title(uniDrugEvo{g},'Color',ABXmech.Colors{tInd},'FontName','Avenir Black','FontSize',16)
    %legend('Line 1', 'Line 2', 'Line 3', 'Line 4', 'Line 5', 'Location','southoutside', 'NumColumns', 2)
    spI = spI+1;  
end

exportgraphics(fig,'Fig2C_ABXRadar.eps','ContentType', 'vector');
ABXorder = uniDrugEvo; 
%% Figure 2B: non-ABX evo Strains: Radar plots per evoDrug: categorized by ABX mechanism 
% get # of cross-resistant ABX for each category 
nonABXMLT = maskLFCtable(maskLFCtable.Type==2,:);
[uniDrugEvo, ia]= unique(nonABXMLT.Drug);
nonABXmechs = [nonABXMLT.Mechanism{ia}];
uniDrugEvo(strcmp(uniDrugEvo, 'Lamotrigine'))= []; % Removing Lamotrigine
mechCrossRes = zeros(length(uniDrugEvo),8); % getting # of cross resistance per mechanism 
mechNumInds={};
for k = 1:length(uniMech)
     for g = 1:length(uniDrugEvo)
             mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),mechInds{k}};
             mechCrossRes(g,k) = sum((mydata>0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
     end
end

mechNames= {"50s","30s","nitrofuran","DNA","folic acid","RNA",sprintf('cell      \n membrane'),"cell wall"}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);


fig = figure('color','w','Position', get(0, 'Screensize')); 
        spI = 1; 
        r=4; NumEvoLines=6;
for g = 1:length(uniDrugEvo)
    data = mechCrossRes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    subplot(5,5,spI)
    
    % spider plot
    spider_plot(data,  'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 6,'LabelFont', 'Avenir Black','LabelFontSize', 8,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,uniDrugEvo{g}));
    title(uniDrugEvo{g},'FontName','Avenir Black','FontSize',14)
    %legend('Line 1', 'Line 2', 'Line 3', 'Line 4', 'Line 5', 'Location','southoutside', 'NumColumns', 2)
    spI = spI+1;  
end
exportgraphics(fig,'Fig2C_nonABXRadar.eps','ContentType', 'vector');
nonABXorder = uniDrugEvo; 
save("../FigureOrder.mat", "nonABXorder", "ABXorder")
%% Supplementary Figure 3:Controls
% get # of cross-resistant ABX for each category 
uniDrugEvo = {"Control", "Water Control"}; 
mechCrossRes = zeros(length(uniDrugEvo),8); % getting # of cross resistance per mechanism 
mechNumInds={};
for k = 1:length(uniMech)
    mechInds{k} = ABXmech.Labels((contains(ABXmech.Mechanism, uniMech{k})));
     for g = 1:length(uniDrugEvo)
             mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),mechInds{k}};
             mechCrossRes(g,k) = sum((mydata>0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
     end
end

mechNames= {sprintf("50s\n(3 drugs)"),sprintf("30s\n(2 drugs)"),sprintf("nitrofuran\n(2 drugs)"),sprintf("DNA\n(3 drugs)"),sprintf("folic acid\n(3 drugs)"),sprintf("RNA\n(2 drugs)"),sprintf('cell membrane\n (3 drugs)'),sprintf("cell wall\n(3 drugs)")}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);

% Resistant 
fig = figure('color','w'); 
tiledlayout(2,2)
spI = 1; 
for g = 1:length(uniDrugEvo)
    nexttile
    data = mechCrossRes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    
    % spider plot
    spider_plot(data, 'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 12,'LabelFont', 'Avenir','LabelFontSize', 10,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,uniDrugEvo{g}));

        title(sprintf("%s", uniDrugEvo{g}),'FontSize',16, 'FontName','Avenir','Color','k')
    spI = spI+1;  
end

% Sensitive 
mechCrossSes = zeros(length(uniDrugEvo),8); % getting # of cross resistance per mechanism 
mechNumInds={};
for k = 1:length(uniMech)
    mechInds{k} = ABXmech.Labels((contains(ABXmech.Mechanism, uniMech{k})));
     for g = 1:length(uniDrugEvo)
             mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),mechInds{k}};
             mechCrossSes(g,k) = sum((mydata<0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
     end
end


for g = 1:length(uniDrugEvo)
    nexttile
    data = mechCrossSes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    
    % spider plot
    spider_plot(data, 'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 12,'LabelFont', 'Avenir','LabelFontSize', 10,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,uniDrugEvo{g}));

        title(sprintf("%s", uniDrugEvo{g}),'FontSize',16, 'FontName','Avenir','Color','k')
    spI = spI+1;  
end
exportgraphics(fig,'Fig2suppControl.svg','ContentType', 'vector');
%% Supplementary Fig 4: Sensitive: ABX evo Strains: Radar plots per evoDrug: categorized by ABX mechanism 
% get # of cross-resistant ABX for each category 
% if want to start 0 in middle change axes offset to 0
abxLFC = maskLFCtable(maskLFCtable.Type==1,:); 

[uniDrugEvo, uniInds] = unique(abxLFC.Drug(abxLFC.Type==1));
mechUni = string(abxLFC.Mechanism(uniInds));
desireMechOrder = ["50s", "30s", "nitrofuran", "DNA", "folicacid","RNA","cellmembrane", "cellwall"];
orgInds = [];
for j9=1:length(desireMechOrder) 
    orgInds = [orgInds; find(contains(mechUni, desireMechOrder(j9)))];
end
uniDrugEvo = uniDrugEvo(orgInds); % REORDER via mechanism in same order as radar plot
mechCrossSes = zeros(length(uniDrugEvo),8); % getting # of cross sensitive permechanism 
for k = 1:length(uniMech)
    for g = 1:length(uniDrugEvo)
        mydata = abxLFC{(contains(abxLFC.Drug,uniDrugEvo{g})),mechInds{k}};
        mechCrossSes(g,k) = sum((mydata<0),'all')/numel(mydata); % getting total # of sensitive drugs in category across lines
    end
end

mechNames= {"50s","30s","nitrofuran","DNA","folic acid","RNA",sprintf('cell      \n membrane'),"cell wall"}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);
fig = figure('color','w','Position', get(0, 'Screensize')); 
        spI = 1; 
        r=4; NumEvoLines=6;
for g = 1:length(uniDrugEvo)
    data = mechCrossSes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    subplot(4,5,spI)
    
    % spider plot
     spider_plot(data,  'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 6,'LabelFont', 'Avenir Black','LabelFontSize', 8,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,erase(uniDrugEvo{g}, " ")));
    title(uniDrugEvo{g},'Color',ABXmech.Colors{tInd},'FontName','Avenir Black','FontSize',16)
    %legend('Line 1', 'Line 2', 'Line 3', 'Line 4', 'Line 5', 'Location','southoutside', 'NumColumns', 2)
    spI = spI+1;  
end

exportgraphics(fig,'Fig2_crossensABX.eps','ContentType', 'vector');
ABXorder = uniDrugEvo; 
%% Supplementary Fig 4: Sensitive: non-ABX evo Strains: Radar plots per evoDrug: categorized by ABX mechanism 
% get # of cross-resistant ABX for each category 
nonABXMLT = maskLFCtable(maskLFCtable.Type==2,:);
[uniDrugEvo, ia]= unique(nonABXMLT.Drug);
nonABXmechs = [nonABXMLT.Mechanism{ia}];
uniDrugEvo(strcmp(uniDrugEvo, 'Lamotrigine'))= []; % Removing Lamotrigine
mechCrossSes = zeros(length(uniDrugEvo),8); % getting # of cross resistance per mechanism 
mechNumInds={};
for k = 1:length(uniMech)
     for g = 1:length(uniDrugEvo)
             mydata = maskLFCtable{(contains(maskLFCtable.Drug,uniDrugEvo{g})),mechInds{k}};
             mechCrossSes(g,k) = sum((mydata<0),'all')/numel(mydata); % getting total # of resistance drugs in category across lines
     end
end

mechNames= {"50s","30s","nitrofuran","DNA","folic acid","RNA",sprintf('cell      \n membrane'),"cell wall"}; 
mechColors=[113 72 157; 158 93 165; 238 126 128; 219 30 62; 12 128 64; 245 138 31; 43 58 150; 252 214 16]./255;

axLim = zeros(2,8);
axLim(2,:) = ones(1,8);


fig = figure('color','w','Position', get(0, 'Screensize')); 
        spI = 1; 
        r=4; NumEvoLines=6;
for g = 1:length(uniDrugEvo)
    data = mechCrossSes(g,:);
    [~, I] = max(data); 
    if sum(data<0.5)==numel(data)
        fColor = [160/255, 160/255, 160/255];
        v=0;
    elseif sum(data>=0.5)>=3
        fColor = [0 0 0];
        v= 0.2;
    else
        fColor = mechColors(I,:);
        v=0.2;
    end 
    subplot(5,5,spI)
    
    % spider plot
    spider_plot(data,  'AxesLabels', 'none', 'AxesDisplay', 'none','AxesInterval',4,'AxesLimits', axLim, ...
        'Color', fColor, 'AxesFontSize', 6,'LabelFont', 'Avenir Black','LabelFontSize', 8,'FillOption', 'on', 'FillTransparency', v, 'AxesPrecision', 2, 'BackgroundColor', 'w', 'AxesLabelsEdge', 'none');
    tInd = find(contains(ABXmech.Labels,uniDrugEvo{g}));
    title(uniDrugEvo{g},'FontName','Avenir Black','FontSize',14)
    %legend('Line 1', 'Line 2', 'Line 3', 'Line 4', 'Line 5', 'Location','southoutside', 'NumColumns', 2)
    spI = spI+1;  
end
exportgraphics(fig,'Fig2_crossensNonABX.eps','ContentType', 'vector');
nonABXorder = uniDrugEvo; 
save("../FigureOrder.mat", "nonABXorder", "ABXorder")