%% Figure 1 & Supplementary Figure 1-2: Evolution vs Screen Paper
%% Load data
cd ../
load("AutoResistanceData_2026-04-24.mat")
load("ABXmechColorMap.mat")
load("allLFC_021926.mat")
load("FigureOrder.mat")
cd Figure1_SupFig1-2/
ICtimepoint = 8;
%% Figure 1C:Bar graph of IC change from wt + Fig 1b: Pie chart of mechanisms 
meanLFCICs = nan(length(uniqueDrugs),1);
ttestICs =  nan(length(uniqueDrugs),1);
meanLineLFCICs = cell(length(uniqueDrugs),1);
types = cell(length(uniqueDrugs),1);
mechs = cell(length(uniqueDrugs),1);
wtICs = nan(length(uniqueDrugs),1);
lineICs = cell(length(uniqueDrugs),1);
for d = 1:length(uniqueDrugs)
    idx = find(IC50table.Timepoint==ICtimepoint & strcmp(uniqueDrugs(d),IC50table.Drug));
    strains = IC50table.FullName(idx);
    wtidx = strcmp(strains, 'MiNoLi wt');
    linesIC = cell2mat(IC50table.RepsIC50(idx(~wtidx)));
    wtICs(d) = IC50table.AvgIC50(idx(wtidx));
    meanLFCICs(d)= mean(linesIC,'omitnan')./wtICs(d); % average LFC of all replicates of each line
    lineICs{d} = IC50table.AvgIC50(idx(~wtidx));
    meanLineLFCICs{d} = IC50table.AvgIC50(idx(~wtidx))./wtICs(d); % average LFC of 3 replicates of each line
    % mean of Lines > mean of wt: using unequal because variances are
    % unlikely to be the same
    ttestICs(d) = ttest2(linesIC,wtICs(d),'Tail', 'right','Alpha',0.05, 'Vartype','unequal');
    types(d) = unique(allData.Type(strcmp(uniqueDrugs(d), allData.DrugName) & strcmp("MiNoLi wt", allData.Strain)));
    mechs(d) = unique(allData.Mechanism(strcmp(uniqueDrugs(d), allData.DrugName) & strcmp("MiNoLi wt", allData.Strain)));
end


fig = figure('color','w');
type1 = ["ABX", "nonABX"];
for i=1:2
    idx= find(strcmp(type1(i),types));
    typemech = mechs(idx);
    x = 1;
    indexes = [];
    logMean = log2(meanLFCICs(idx));
    [sortLogMean, sIdx] = sort(log2(meanLFCICs(idx)), 'descend');
    set(0, 'CurrentFigure', fig)
    subplot(2,1,i)
    hold on;
    temp1 = [];
    for i2 = sIdx'
        mechidx =find(strcmp(rowColors.Labels, typemech(i2)));
        bar(x, logMean(i2), 'FaceColor',rowColors.Colors{mechidx})
        x=x+1;
        temp1 = [temp1, logMean(i2)];
        indexes = [indexes; idx(i2)];
    end

    d2=1;temp2 =[];
    for d1 = indexes'
        scatter(zeros(length(meanLineLFCICs{d1}),1)+d2, log2(meanLineLFCICs{d1}),'ko','filled','XJitter','rand','XJitterWidth',0.5)
        temp2 = [temp2; log2(meanLineLFCICs{d1})'];
        if ttestICs(d1) ==1
            plot(d2, log2(meanLFCICs(d1))+1, "r*")
        end
        d2 = d2+1;
    end
    grid on;box on;
    xticks(1:x)
    ylim([-2 12])
    xticklabels(uniqueDrugs(indexes))
    ylabel('L2FC of IC50 (uM)')
    title(sprintf('%s evolved drugs',type1(i)))
    if i==1
        % Saving ABX LFC infor
        LFCIC = temp1;
        LFCICpoints = temp2; 
        labelIC = uniqueDrugs(indexes);
    end 
end
sgtitle(sprintf("%i hr", ICtimepoint))
saveas(fig, 'Fig1D_IC50bar_8hr.png')
exportgraphics(fig, 'Fig1D_IC50bar_8hr.svg', 'ContentType', 'vector');

%% Supplementary Figure 1: Comparing autoresistance between IC and BMD bar
abxLFC = maskLFCtable(maskLFCtable.Type==1,:); 
labelBMD = unique(abxLFC.Drug); 
labelBMD{strcmp(labelBMD,'Tetracycline HCl')} = 'Tetracycline';
autoRes = nan(length(labelBMD),4);
mechidx = nan(length(labelBMD),1);
for i = 1:length(labelBMD)
    mechidx(i) =find(contains(ABXmech.Labels, labelBMD{i}));
    autoRes(i, :) = abxLFC{contains(abxLFC.Drug, labelBMD{i}), contains(abxLFC.Properties.VariableNames, labelBMD{i})};
end 
LFCBMD = mean(autoRes,2);

% Comparing autoresistance between IC and BMD as scatter
figure 
hold on
storeLFCBMD = [];
for p1 = 1:length(labelIC)
    if strcmp(labelIC(p1),'Tetracycline HCl')
       idx = find(contains(labelBMD, 'Tetracycline'));
    else 
        idx = find(contains(labelBMD, labelIC(p1)));
    end 
    scatter(LFCIC(p1),LFCBMD(idx), 'filled','MarkerFaceColor', ABXmech.Colors{mechidx(idx)})
    text(LFCIC(p1)+0.025,LFCBMD(idx), labelIC(p1), 'HorizontalAlignment','left')
    storeLFCBMD = [storeLFCBMD, LFCBMD(idx)];
end 
grid on; box on;
ylim([0,8])
xlabel('log2FC(IC)')
ylabel('log2FC(BMD)')
p = polyfit(LFCIC, storeLFCBMD,1);
plot(0:8, polyval(p,0:8), 'r')
[R,Pvals] = corrcoef(LFCIC, storeLFCBMD);
title(sprintf('r = %f pvalue = %f', R(1,2), Pvals(1,2)))
plot(0:8,0:8, 'k--')

% Comparing autoresistance between IC and BMD bar
fig2 = figure;
subplot(2,1,1)
hold on
storeLFCBMD = [];
for p1 = 1:length(labelIC)
    if strcmp(labelIC(p1),'Tetracycline HCl')
       idx = find(contains(labelBMD, 'Tetracycline'));
    else 
        idx = find(contains(labelBMD, labelIC(p1)));
    end 
    bar(p1,LFCIC(p1),'FaceColor', ABXmech.Colors{mechidx(idx)})
    scatter(zeros(1,4)+p1, LFCICpoints(p1,:),'ko','XJitter','rand','XJitterWidth',0.5)
end 
xticks(1:length(labelIC))
xticklabels(labelIC)
grid on; box on;
ylim([0,8])
ylabel('log2FC')
title(sprintf('Pos=IC; Neg=BMD; r = %f pvalue = %f', R(1,2), Pvals(1,2)))

subplot(2,1,2)
hold on
storeLFCBMD = [];
for p1 = 1:length(labelIC)
    if strcmp(labelIC(p1),'Tetracycline HCl')
       idx = find(contains(labelBMD, 'Tetracycline'));
    else 
        idx = find(contains(labelBMD, labelIC(p1)));
    end 
    bar(p1, LFCBMD(idx), 'FaceColor', ABXmech.Colors{mechidx(idx)})
    scatter(zeros(1,4)+p1, autoRes(idx,:),'ko','XJitter','rand','XJitterWidth',0.5)
end 
xticks(1:length(labelIC))
xticklabels(labelIC)
grid on; box on;
ylim([0,8])
ylabel('log2FC')
title('BMD')
saveas(fig2, 'Fig1Supplementary_ICvsBMD.png')
exportgraphics(fig2, 'SupplementaryFig2_ICvsBMD.svg', 'ContentType', 'vector');
