%% AutoResistance_ICdata: To calculate IC50 from dose-response data 
load("Evo_IC50Data_2026-04-22.mat")
ICpoint = 50; 
ICtimepoints = [8];
timeInterval = 20/60;
%% Plot Median of Normalize Growth curves as sanity checks
% Using 'blankedData'= subtract the average of 0-40mins for each well 
% Plotting by strain
qSP2=1;counter = 1;
uniqueDrug = unique(allData.DrugName);
for dn1 = 1:length(uniqueDrug)
    if rem(qSP2,25) == 1
        qSP2 = 1;
    fig4 = figure('color','w','Position', get(0, 'Screensize'));
    end 
    drugIDX = strcmp(uniqueDrug(dn1), allData.DrugName);
    uniqueDilutions = unique(allData.uM(drugIDX))';
    c = [sort(uniqueDilutions,'descend')',turbo(length(uniqueDilutions))];
    uniwt = unique(allData{drugIDX,"Strain"});
    uniwt = uniwt([5 1:4]);
    for l =1:length(uniwt)
        set(0,'CurrentFigure',fig4)
        subplot(5,5,qSP2)
        qSP2 = qSP2 + 1;
        hold on
        DataIndex = find(strcmp(uniwt(l),allData.Strain) & strcmp(uniqueDrug(dn1), allData.DrugName));
        [~, inds] =sort(allData.uM(DataIndex));
        sortDataIndex = DataIndex(inds);
        
        plotsVector = [];
        dilutions = [];
        % Plotting each strain
        for kp = 1:length(DataIndex)
            md = allData.BlankedData{sortDataIndex(kp)};
            time = 0:timeInterval:(length(md)-1)*timeInterval;
            dil = allData.uM(sortDataIndex(kp));
            myind = find(dil==c(:,1));
            for kp1 = 1:width(md)
                if kp1==1
                p= plot(time, md(:,kp1),'Color', c(myind,2:end),'LineWidth', 2);
                else
                    plot(time, md(:,kp1),'Color', c(myind,2:end), 'LineWidth', 2);
                end 
            end 
            if ~ismember(dil,dilutions)
                plotsVector = [plotsVector, p];
                dilutions = [dilutions , dil];
            end
        end
        xline(8)
        xlabel('Time(hr)')
        ylabel('Abs(OD600)')
        ylim([-.1 1.2])
        title(sprintf("%s %s", uniwt{l}, uniqueDrug{dn1}))
        legend(plotsVector, string(dilutions), 'Location', 'eastoutside')
        hold off
        grid on
    end
    if rem(qSP2,25) == 1
        sgtitle(sprintf("Part %i Normalized Growth Curves",counter))
        saveas(fig4,sprintf('./NormalizedGrowthCurves/NormalizedGrowthCurveConc_Part%i.jpg',counter ))
        counter = counter+1;
    end 
end

%% Plot IC50 curves 
fprintf('Finding IC50 \n')
plotNum = 1;
% Concatenating all plate data together per experimental run
% Control wells in all plates for experimental run will be averaged together
dilution = allData.uM;
DrugType = allData.DrugName;
strains = allData.Strain;
uniqueDrugs = unique(DrugType);
rowIndex = 1;
IC50table = {};
strainNum = 1;

time = 0:timeInterval:72*timeInterval;
counter = 0; r= 3; c =6;
for dn2 = 1:length(uniqueDrugs)
    for tp1 = 1:length(ICtimepoints)
        if rem(plotNum, r*c) == 1
            if counter>0
            sgtitle(sprintf('Set%i: IC%i Growth Curves', counter, ICpoint))
            saveas(fig2,sprintf('./IC50/IC%i_Curves_Set%i.jpg', counter, ICpoint))
            end
            counter = counter + 1;
            fig2 = figure('color','w','Position',  get(0, 'Screensize'));
            plotNum = 1;
        end 
        subplot(r,c,plotNum)
        plotNum = plotNum + 1;
        if strcmp(uniqueDrugs(dn2), "PBS") % Skipping PBS wells
            continue
        end
        strainTypes = unique(strains(strcmp(uniqueDrugs(dn2), DrugType)));
        strainTypes(strcmp(strainTypes,'PBS'))=[];
        hold on
        tpoints = [];
        legendPlots = [];
        colors = parula(length(strainTypes));
        for strain = 1:length(strainTypes)
            IC50array = [];
            graphCurves = {};
            predCurves = {};
            coefficients = {};
            ind = find(strcmp(strainTypes{strain}, strains) & strcmp(uniqueDrugs(dn2), DrugType));
            if isempty(ind)
                continue %Skipping drugs that was not tested with the strain
            end
            strainNum = strainNum+1; % adding to strain number here so it won't be added to when "continue" is run
            dil = dilution(ind);
            normData = allData.BlankedData(ind);
            timeIDX = find(time==(ICtimepoints(tp1)));
            maxGrowth = cellfun(@(x) x(timeIDX,:), normData, 'UniformOutput', false);

            % deal with 0 dosage by using it to normalise the results
            response = [];
            dose = [];
            controlResponse = mean(maxGrowth{dil == 0},"omitnan"); % average of no drug
            for mi = 1:length(dil)
                if dil == 0
                    dilresponse = ones(1,length(maxGrowth{dil == 0}));
                else
                    dilresponse = maxGrowth{mi}/controlResponse;
                    numReps = 3;
                    if length(dilresponse)<numReps % Adding in nan for when replicates are missing
                        if isscalar(dilresponse)
                            x= NaN(numReps,1);
                            x(1) = dilresponse;
                            dilresponse = x;
                        else
                            dilresponse = padarray(dilresponse, [0 numReps-length(dilresponse)], nan, 'post');
                        end
                    end
                end
                dose = [dose; dil(mi)];
                response =[response; dilresponse];
            end

            ypointsPlot = [];
            for replicates = 1:width(response) % to handle situtations where replicates is more than 3
                currentResponse = response(:,replicates); % Response of the replicate being analyzed
                % Using aICarray for absolute IC50 because for some cases
                % did not reach a lower plateau 
                [coeffs,rICarray,aICarray,H1,H2, xpoints, ypoints] = calcDoseResponse_v2(dose,currentResponse,false ,ICpoint);
                IC50 = aICarray(2);
                responseAtIC = aICarray(3);
                % Handling the situations where the strain is very resistant
                % and the response is greater then without drug so IC50 is
                % off the chart but will be made to be maximum dose tested 
                if  currentResponse(max(dose)==dose)>1 % If responses at the highest dose is greater than 1
                    IC50 = max(dose);
                elseif IC50 <= 0 % When IC50<=0 because the strain is very resistant
                    IC50 = max(dose);
                elseif currentResponse(length(currentResponse)-1)>1 && IC50 < dose(length(currentResponse)-1)
                    % When responses at 2nd highest dose is greater than 1
                    % and IC50< 2nd highest dose (this would be because the
                    % IC50 is taken from a positive slope)
                    IC50 = max(dose);
                elseif isnan(IC50)
                    IC50 = max(dose);
                end
                IC50array(replicates,:) = [IC50, responseAtIC];
                graphCurves{replicates} = {dose,currentResponse};
                predCurves{replicates} = {xpoints,ypoints};
                ypointsPlot = [ypointsPlot, ypoints'];
                coefficients{replicates} = coeffs;
            end

            % Graphing IC50 curves
            c1 = colors(strain,:);
            p3 = plot(IC50array(:,1),IC50array(:,2), "*",'LineWidth', 4, 'Color', c1);
            legendPlots=[legendPlots,p3];
            p =plot(dose, response,  'o','Color',c1);
            p2 = plot(xpoints, ypointsPlot, 'LineWidth', 2, 'Color',c1);


            title(sprintf("%s: %i hr", uniqueDrugs{dn2}, ICtimepoints(tp1)))
            xlabel('Dose (uM)')
            ylabel('normalized abosrbance 600nm')
            ylim([-0.1,1.5])
            xlim([0, max(dose)])
            set(gca, 'Xscale', 'log')
            tpoints = [tpoints , sprintf("%i H %s IC%i: %0.2f", ICtimepoints(tp1), strainTypes{strain}, ICpoint, mean(IC50array(:,1), 'omitnan'))];

            IC50table{rowIndex,1} = uniqueDrugs(dn2);
            IC50table{rowIndex,2} = IC50array(:,1); % relative IC50
            IC50table{rowIndex,3} = mean(IC50array(:,1), 'omitnan');
            IC50table{rowIndex,4} = graphCurves'; %Storing dose and response
            IC50table{rowIndex,5} = strainTypes{strain}; %Storing full name
            IC50table{rowIndex,6} = max(dil); %storing max conc
            IC50table{rowIndex,7} = predCurves';
            IC50table{rowIndex,8} = coefficients; %storing coefficients for sigmoid equation
            IC50table{rowIndex,9} = ICtimepoints(tp1); %storing timepoint
            IC50table{rowIndex,10} = IC50array(:,2); % response to relative IC50
            rowIndex = rowIndex +1;
        end
        grid on
        legend(legendPlots,tpoints, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 1 )
    end
end
%Saving IC50 Information
IC50table = cell2table(IC50table);
IC50table.Properties.VariableNames{'IC50table2'} = 'RepsIC50'; 
IC50table.Properties.VariableNames{'IC50table3'} = 'AvgIC50';
IC50table.Properties.VariableNames{'IC50table4'} = 'DoseResponse';
IC50table.Properties.VariableNames{'IC50table5'} = 'FullName';
IC50table.Properties.VariableNames{'IC50table6'} = 'MaxConc';
IC50table.Properties.VariableNames{'IC50table7'} = 'PredictedCurve';
IC50table.Properties.VariableNames{'IC50table8'} = 'Coefficients';
IC50table.Properties.VariableNames{'IC50table1'} = 'Drug';
IC50table.Properties.VariableNames{'IC50table9'} = 'Timepoint';
IC50table.Properties.VariableNames{'IC50table10'} = 'ReponseRepsIC50';

fprintf('Finished finding IC50 \n')
%% 
save("AutoResistanceData_2026-04-24.mat",'IC50table','uniqueDrugs', 'allData') 
