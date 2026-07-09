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
ttestICs =  nan(length(uniqueDrugs),2);
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
    [ttestICs(d,1), ttestICs(d,2)] = ttest2(linesIC,IC50table.RepsIC50{idx(wtidx)},'Tail', 'right','Alpha',0.05, 'Vartype','unequal');
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
        if ttestICs(d1,1)==1 && ttestICs(d1,2)<0.0001
            plot(d2, log2(meanLFCICs(d1))+1.5, "k*")
        elseif ttestICs(d1,1)==1 && ttestICs(d1,2)<0.001
            plot(d2, log2(meanLFCICs(d1))+1.5, "g*")
        elseif ttestICs(d1,1)==1 && ttestICs(d1,2)<0.01
            plot(d2, log2(meanLFCICs(d1))+1, "b*")
        elseif ttestICs(d1,1)==1 && ttestICs(d1,2)<0.05
            plot(d2, log2(meanLFCICs(d1))+1, "r*")
        else
            plot(d2, log2(meanLFCICs(d1))+1, "ro")
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
        % Saving ABX LFC info
        LFCIC = temp1;
        LFCICpoints = temp2; 
        labelIC = uniqueDrugs(indexes);
    else
        totLFCIC = [LFCIC, temp1];
        totLFCICpoints = [LFCICpoints; temp2]; 
        totlabelIC = [labelIC; uniqueDrugs(indexes)];
    end 

    p=[plot(nan,nan, 'k*'),plot(nan,nan, 'g*'),plot(nan,nan, 'b*'),plot(nan,nan, 'r*')];
    legend(p, {'p<0.0001','p<0.001', 'p<0.01','p<0.05'})
end
sgtitle(sprintf("%i hr", ICtimepoint))
saveas(fig, 'Fig1D_IC50bar_8hr.png')
exportgraphics(fig, 'Fig1D_IC50bar_8hr.svg', 'ContentType', 'vector');
%% Saving variables
autoResTable = table(totlabelIC, totLFCIC', totLFCICpoints, 'VariableNames',["drug", "avgLFC_IC", "strainAvgLFC_IFC"]);
cd ../
save("AutoResistanceData_2026-04-24.mat", "autoResTable", '-append')
cd Figure1_SupFig1-2/
%% Supplementary Figure 1: Comparing autoresistance between IC and BMD scatter
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
% LFC to FC
FCICpoints = LFCICpoints;
orderFCIC= [];
FCBMDpoints = autoRes; 
orderFCBMD = [];
for p1 = 1:length(labelIC)
    if strcmp(labelIC(p1),'Tetracycline HCl')
       idx = find(contains(labelBMD, 'Tetracycline'));
    else 
        idx = find(contains(labelBMD, labelIC(p1)));
    end 
    for p2=1:4
        scatter(FCBMDpoints(idx,:),FCICpoints(p1,:),[], ABXmech.Colors{mechidx(idx)},'filled')
    end 
    orderFCBMD = [orderFCBMD, FCBMDpoints(idx,:)];
    orderFCIC= [orderFCIC,FCICpoints(p1,:) ];
    %text(FCICpoints(p1,:),FCBMDpoints(idx,:), labelIC(p1))
end 


yticks(0:8); yticklabels(["1", "2", "4", "8", "16", "32", "64", "128", "256"])
xticks(0:5); xticklabels(["1", "2", "4", "8", "16", "32"])
grid on; box on;
xlim([0,5]), ylim([0,8])
ylabel('IC50 fold change')
xlabel('MGC fold change')
[p,S] = polyfit(orderFCBMD, orderFCIC,1);
plot(0:140, polyval(p,0:140), 'k')
[R,Pvals] = corrcoef(orderFCBMD,orderFCIC);
title(sprintf('r = %f pvalue = %s', R(1,2), Pvals(1,2)))
axis square;



%% Supplementary Figure 2: t-test comparision of auto-resistance by mechanisms
figure; 
tiledlayout(2,1)
uniMech = {'50s', '30s', 'cell membrane', 'cell wall', 'DNA', 'folic acid', 'nitrofuran', 'RNA'};
nexttile 
hold on;
x=[]; g= [];
allSamp = reshape(totLFCICpoints,1,[]);
% Plotting each ABX mechanism 
for u1 = 1:length(uniMech)
    idx = find(strcmp(allData.Mechanism, uniMech{u1}));
    drugsInterest = unique(allData{idx, "DrugName"});
    sbIDX = find(contains(totlabelIC, drugsInterest));
    notIDX = 1:length(totLFCIC);
    notIDX(sbIDX) = [];
    samp = reshape(totLFCICpoints(sbIDX,:), 1,[]); 
    x= [x,samp];
    g = [g;repmat({sprintf('%s(%i)',uniMech{u1}, numel(drugsInterest))},length(samp),1)];
    scatter(zeros(length(samp),1)+u1, samp,'.k', 'XJitter','rand', 'XJitterWidth',0.15)
    yline(median(allSamp), '--')
    [~,p]=ttest(samp,median(allSamp));
    if p<0.0001
        plot(u1,max(samp)+1, 'k*')
    elseif p<0.001
        plot(u1,max(samp)+1, 'g*')
    elseif p<0.01
        plot(u1,max(samp)+1, 'b*')
    elseif p<0.05
        plot(u1,max(samp)+1, 'r*')
    end 
end 
% Plotting nonABX 
uniNonABX = unique(allData{strcmp(allData.Type, 'nonABX'), 'DrugName'});
nonABXidx = find(contains(totlabelIC,uniNonABX));
nonABXvalues = reshape(totLFCICpoints(nonABXidx,:), 1,[]);
[~,p]=ttest(nonABXvalues,median(allSamp));
if p<0.0001
    plot(u1+1,max(nonABXvalues)+1, 'k*')
elseif p<0.001
    plot(u1+1,max(nonABXvalues)+1, 'g*')
elseif p<0.01
    plot(u1+1,max(nonABXvalues)+1, 'b*')
elseif p<0.05
    plot(u1+1,max(nonABXvalues)+1, 'r*')
end
g = [g;repmat({sprintf('nonABX(%i)', numel(uniNonABX))},length(nonABXvalues),1)];
x = [x, nonABXvalues];
scatter(ones(length(nonABXvalues),1)+u1, nonABXvalues,'.k', 'XJitter','rand', 'XJitterWidth',0.4)
% Plotting all 
g = [g;repmat({'all'},length(allSamp),1)];
x = [x, allSamp];
scatter(ones(length(allSamp),1)+u1+1, allSamp,'.k', 'XJitter','rand', 'XJitterWidth',0.4)

% Plotting box plot 
boxplot(x,g)

title("ttest to mean of mechanism to median of all data")
grid on;
xlim([0, u1+3])
p1=[plot(nan,nan, 'k*'),plot(nan,nan, 'g*'),plot(nan,nan, 'b*'),plot(nan,nan, 'r*')];
legend(p1, {'p<0.0001','p<0.001', 'p<0.01','p<0.05'})
ylabel('IC50 change log2 fold change')

% Supplementary Figure: t-test comparision of auto-resistance for multidrug resistance (MDR)
superBugs = {'Chloramphenicol', 'Tetracycline HCl', 'Ciprofloxacin','Mitomycin C', 'Cisplatin','Pentamidine','Dichlorophene','Sertraline'};
sbIDX = nan(length(superBugs),1);
for l1 = 1:length(superBugs)
    temp = find(contains(totlabelIC,superBugs{l1}));
    if ~isempty(temp)
        sbIDX(l1)=temp;
    end
end 

nexttile; hold on; 
notIDX = 1:height(totLFCICpoints);
notIDX(sbIDX) = [];
samp = reshape(totLFCICpoints(sbIDX,:), 1,[]); 
allSamp = reshape(totLFCICpoints(notIDX,:),1,[]);
x= [samp,allSamp];
g = [repmat({sprintf('MDR(%i)', numel(superBugs))},length(samp),1);repmat({'non'},length(allSamp),1)];
[h,p]=ttest2(samp,allSamp);
scatter(ones(length(samp),1), samp,'k.', 'XJitter','rand', 'XJitterWidth',0.15)
scatter(ones(length(allSamp),1)+1, allSamp,'k.', 'XJitter','rand', 'XJitterWidth',0.15)
boxplot(x,g)
ylabel('IC50 change log2 fold change')
title(sprintf("2 sample ttest: h= %i, p= %f", h,p))
grid on; 