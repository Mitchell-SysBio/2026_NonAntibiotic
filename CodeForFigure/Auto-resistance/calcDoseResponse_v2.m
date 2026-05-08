function [coeffs,relICarray,absICarray,fh1,fh2, xpoints, ypoints] = calcDoseResponse_v2(drugConc,normAbsorbance,tfPlot,requestedICs)
% function [coeffs,relICarray,absICarray,fh1,fh2, xpoints, ypoints] = calcDoseResponse(drugConc,normAbsorbance,tfPlot,requestedICs)

% doseresponse function
% Use this function to calculate relative and absolute IC/EC values: 
% Computes the Hill Coefficient and IC50 of a
% dose/response relationship given two vectors describing the doses and
% responses. A semilog graph is plotted illustrating the relationship. The IC50 is also
% labelled. Requires nlinfit from statistics toolbox.
% IC: inhibitory concentrations 
% EC: effective concentrations
% Relative: concentration required to bring out the value you
% want between the top and bottom plateaus of the curve 
% Absolute:concentration required to bring out the value you want between
% the 0 and maximum response 

% Inputs: UPDATE
% drugConc: Provide an array of drug concentrations used. Must be horizontal. If you include a
% zero concentration, it will be replaced by four folds below.
% norm_abs: Normalized(Blank Subtracted and divided by no drug) Absorbance 600nm values. Must be horizontal.  
% plotOpt: If true a plot will be returned
% ICArrayRequested: Enter a vector of desired IC values (e.g. IC25, IC75 -> [25,75])

% Outputs: 
% coeffs: Coefficients of the fit
% rICarray: Relative IC Array
%   - First row: Desired IC value e.g. IC25->25
%   - Second row: IC value concentration
%   - Third row: response at IC value 
% aICarray: Absolute IC Array
%   - First row: Desired IC value e.g. IC25->25
%   - Second row: IC value concentration
%   - Third row: response at IC value 
% H1 and H2: figure handles 
% xpoints, ypoints: values used to plot the sigmoid curve

% Note Version 2: 2024-09-10: Edited by Carmen to round to 4 decimal places for
% drugs that have very low concentrations 
%% Fitting the sigmoid equation

% checking if drugConc and norm_abs is vertical and transposing if it is
% not
[r, ~] = size(drugConc);
if r ~=1
    drugConc = drugConc'; 
end 
[r1, ~] = size(normAbsorbance);
if r1 ~=1
    normAbsorbance = normAbsorbance'; 
end 

% hill equation sigmoid (will be used for the dose curve fitting)
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

% dealing with drug conc 0 (0 conc. will be replaced with very small conc.)
if ~isempty(find(drugConc==0,1))
   drugConcSorted = sort(drugConc,'ascend');
   fc = drugConcSorted(1,3)./drugConcSorted(1,2); % infer the conc jump by 2nd and 3rd conc.
   drugConc(1,drugConc==0) = drugConcSorted(1,2)./(fc^4); % swap with conc. smaller by 4-fold jumps 
end  

% calculate some rough guesses for initial parameters
minResponse = min(normAbsorbance);
maxResponse = max(normAbsorbance);
midResponse = mean([minResponse maxResponse]);
minConc = min(drugConc);
maxConc= max(drugConc);

% Turning off warning messages
warning('off','all') 

% fit the curve and to a sigmoid function 
[coeffs,r,J]=nlinfit(drugConc,normAbsorbance,sigmoid,[minResponse maxResponse midResponse 0]);
IC50=coeffs(3);
hillCoeff=coeffs(4);

%% calculate requested ICs
NumberXpoints = 5000; % number of points to generate for xpoints 
xpoints=logspace(log10(minConc),log10(maxConc),NumberXpoints); 
ypoints = sigmoid(coeffs,xpoints);

for i=1:length(requestedICs)
        curIC = requestedICs(i);
        relResponseSpan = max(ypoints)-min(ypoints);
        % Calculate relative ICs    
        relResponseAtIC =  ((100-curIC)/100)*relResponseSpan+min(ypoints); 
        [~,relInx] = min(abs(ypoints-relResponseAtIC));
        relCurDose = round(xpoints(relInx),4); % edited 2024-09-10 (Carmen) to allow for more decimal points for really low drug conc
        if relInx == NumberXpoints
            relCurDose = nan;
            relResponseAtIC = nan;
        end 
        relICarray(1,i) = curIC; %requested IC value (ex. 50, 75, etc.)
        relICarray(2,i) = relCurDose; % IC 
        relICarray(3,i) = relResponseAtIC; % response at IC dosage 

        % Calculate absolute ICs
        absResponseSpan = max(ypoints)-0;
        absResponseAtIC =  ((100-curIC)/100)*absResponseSpan+0;
        [~,absInx] = min(abs(ypoints-absResponseAtIC));
        absCurDose = round(xpoints(absInx),2);
        if absInx == NumberXpoints
            absCurDose = nan;
            absResponseAtIC = nan;
        end 
        absICarray(1,i) = curIC;
        absICarray(2,i) = absCurDose;
        absICarray(3,i) = absResponseAtIC;
        
end

%% plot the fitted sigmoid

if(tfPlot)
    %anotate the relative IC requested by user
    fh1=figure('color','white');
    hold on;
    grid on; 
    semilogx(xpoints,ypoints,'Color',[0 0 0],'LineWidth',2);
    ylim([0 1.2]);
    plot(drugConc, normAbsorbance, "b*", "LineWidth", 2);
    xlabel("log(Dose) (uM)");
    ylabel("Normalized Absorbance A600nm");
    text(IC50,mean([coeffs(1) coeffs(2)]),[sprintf('rIC_{50}=%0.2f',IC50) ' \rightarrow '],'FontSize',16,'color',[0.5 0.5 0.5],'HorizontalAlignment','right' );
    title("Dose Reponse: Relative IC50");
    for i=1:length(requestedICs)
        curSTR = ['relIC' num2str(relICarray(1,i)) '=' num2str(relICarray(2,i))];
        text(relICarray(2,i),relICarray(3,i),[' \leftarrow ' sprintf('%s',curSTR)],'FontSize',12);
    end
    %anotate the absolute IC requested by user
    fh2=figure;
    hold on;
    grid on; 
    semilogx(xpoints,ypoints,'Color',[0 0 0],'LineWidth',2);
    plot(drugConc, normAbsorbance, "b*", "LineWidth", 2)
    xlabel("log(Dose) (uM)")
    ylabel("Normalized Absorbance A600nm")
    ylim([0 1.2])
    title("Dose Reponse: Absolute IC50");
    for i=1:length(requestedICs)
       curSTR = ['absIC' num2str(absICarray(1,i)) '=' num2str(absICarray(2,i))];
        text(absICarray(2,i),absICarray(3,i),[' \leftarrow ' sprintf('%s',curSTR)],'FontSize',12);
    end
else
    fh1 = 0;
    fh2 = 0; 
end 

end
