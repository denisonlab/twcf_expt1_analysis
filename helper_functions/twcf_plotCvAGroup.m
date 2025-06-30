function [aGroup] = twcf_plotCvA(varA,varC,p,a)
% [a] = twcf_plotCvA(varA,varC,p,a)
% Plots type-2 function of Consciousness (C) vs Accuracy (A) 
% using type-1 fitted params to C(x) and A(x), to test PSI 
% Inputs:
    % varA: accuracy structure 
    % varC: consciousness structure 
    % p: analysis params (from twcf_analysisParams.m) 
    % a: analysis structure (from behav_analysis.m) 
% Ouputs: 
    % a: analysis structure with PSI analysis 

%% Set up subplot 
subplot 111
subplot('Position', [0.1570 0.2074 0.7176 0.7176])
figureStyle

hold on
axis square

xlim([0.441 1]) % centers the ticks of the chance lines
ylim([0 1])

xticks([0.5 0.75 1])
yticks([0 0.25 0.5 0.75 1])

% chance lines
vl = plot([0.5 0.5],[0 1],'k','LineStyle','--'); % vertical line 
hl = plot([0 1],[0.5 0.5],'k','LineStyle','--'); % horizontal line
uistack(vl,"bottom")
uistack(hl,"bottom")

xlabel('p(correct discrimination)')
ylabel(sprintf('p("%s")',varC.varShortName))

xAX = get(gca,'XAxis');  
yAX = get(gca,'YAxis');  
xlab = get(gca,'XLabel');
ylab = get(gca,'YLabel');

set(xAX,'FontSize', p.style.textTickSize)
set(yAX,'FontSize', p.style.textTickSize)
set(xlab,'FontSize', p.style.textAxisSize)
set(ylab,'FontSize', p.style.textAxisSize)

%% Plot CvA by attention 
for iAtt = 1:numel(varC.attHeaders) 
    % Plot datapoints
    startIdx = 1; 
    switch a(1).exptShortName
        case 'twcf_cue_gab_det'
            for i = 1:numel(varC.stimLevels)
                AX(i) = varA.var(i,iAtt); AX(1)=0.5; 
                CX(i) = varC.var(i,iAtt);
            end
        otherwise
            for i = startIdx:numel(varC.stimLevels)
                AX(i) = varA.var(i,iAtt);
                CX(i) = varC.var(i,iAtt);
            end
    end
    sPlot(iAtt) = scatter(AX,CX,p.style.sz,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',1,'MarkerEdgeColor','w','LineWidth',2);

    %% Plot fits
    clear A; clear C;
    switch a(1).exptShortName
        case 'twcf_cue_gab_det'
            startIdx = 2; % 2; don't plot 0 contrast; nevermind we want this
        otherwise
            startIdx = 1;
    end
    sIdx = varA.stimLevelsFine >= varA.stimLevels(startIdx); 
    A = varA.fit(iAtt,sIdx);
    C = varC.fit(iAtt,sIdx);
    plot(A,C,'Color',p.style.attColors(iAtt,:),'LineWidth',2);
    
    % Constrain CvA by range of observed A and/or C data values 
    % across stimulus levels and attention condition 
    % Otherwise uses theoretical bounds 
    boundA = 1; boundC = 0; 
    if boundC % data-driven bounds
        [rangeC(iAtt,1),rangeC(iAtt,2)] = bounds(CX,'omitnan');
    else % theoretical bounds
        rangeC(iAtt,1) = 0;
        rangeC(iAtt,2) = 1;
    end
    if boundA % data-driven bounds
        % [rangeA(iAtt,1),rangeA(iAtt,2)] = bounds(AX,'omitnan');
        rangeA(iAtt,1) = 0.5; % constrain lower bound regardless 
        % max of fits 
        % rangeA(iAtt,2) = max(varA.fit(iAtt,:)); 
        % max of data 
        rangeA(iAtt,2) = max(varA.var(:,iAtt)); 
        % rangeA(iAtt,2) = max(AX); % end or max? 
    else % theoretical bounds
        rangeA(iAtt,1) = 0.5;
        rangeA(iAtt,2) = 1;
    end 
end

%% Find overlap between C and/or A data across attention conditions
minC = max(rangeC(:,1));
maxC = min(rangeC(:,2));
minA = max(rangeA(:,1));
maxA = min(rangeA(:,2));

% Draw patch showing overlap 
% lc1 = yline(minC,'Color',p.style.labelColors); 
% lc2 = yline(maxC,'Color',p.style.labelColors); 
% la1 = xline(minA,'Color',p.style.labelColors); 
% la2 = xline(maxA,'Color',p.style.labelColors); 

% uistack(lc1,'bottom')
% uistack(lc2,'bottom')
% uistack(la2,'bottom')

clear x; clear y; 
xBound = [minA minA maxA maxA];
yBound = [minC maxC maxC minC];
% patch(xBound,yBound,p.style.labelColors); 
% set(gca,'children',flipud(get(gca,'children'))) % send patch back 

%% Get fitted params by attention 
for iAtt = 1:numel(varC.attHeaders) 
    % C params
    alphaC   = varC.paramsFitted(iAtt,find(contains(varC.paramValueNames,'threshold')));
    betaC    = varC.paramsFitted(iAtt,find(contains(varC.paramValueNames,'slope')));
    gammaC   = varC.paramsFitted(iAtt,find(contains(varC.paramValueNames,'guess-rate')));
    lambdaC  = varC.paramsFitted(iAtt,find(contains(varC.paramValueNames,'lapse-rate')));

    % A params
    alphaA   = varA.paramsFitted(iAtt,find(contains(varA.paramValueNames,'threshold')));
    betaA    = varA.paramsFitted(iAtt,find(contains(varA.paramValueNames,'slope')));
    gammaA   = varA.paramsFitted(iAtt,find(contains(varA.paramValueNames,'guess-rate')));
    lambdaA  = varA.paramsFitted(iAtt,find(contains(varA.paramValueNames,'lapse-rate')));
    
    % Relative params
    rA = alphaC/alphaA; % relative threshold
    rB = betaC/betaA; % relative slope 
    rG = gammaC/gammaA; % relative guess-rate 
    rL = lambdaC/lambdaA; % relative lapse-rate 
    
    % Save 
    % a.alphaC(iAtt) = alphaC; 
    % a.betaC(iAtt) = betaC; 
    % a.gammaC(iAtt) = gammaC; 
    % a.lambdaC(iAtt) = lambdaC; 
    % 
    % a.alphaA(iAtt) = alphaA; 
    % a.betaA(iAtt) = betaA; 
    % a.gammaA(iAtt) = gammaA; 
    % a.lambdaA(iAtt) = lambdaA; 
    % 
    % % CvA param ratios 
    % a.alphaR(iAtt) = rA; 
    % a.betaR(iAtt) = rB; 
    % a.gammaR(iAtt) = rG;
    % a.lambdaR(iAtt) = rL; 

    %% Weibull functions for C and A
    ppd = 0.0202; % edit to be a dynamic input, or 0.07 based on data?, should be 0 stimulus strength (1 px) 
    switch p.fit.PFtype
        case 'gumbel'
            x = log10(ppd:0.01:1.2);
        case 'weibull'
            x = ppd:0.01:1.2;
    end
    Cfit = twcf_fittedPF(alphaC, betaC, gammaC, lambdaC, x, p); 
    Afit = twcf_fittedPF(alphaA, betaA, gammaA, lambdaA, x, p); 
    % plot(Afit,Cfit,':','Color',p.style.attColors(iAtt,:),'LineWidth',2); 

    %% Functional form (Interpreted from Brian's metaperceptual equation) 
    % Uses parameters fitted to stimuli in linear space 
    plotCvAEquation = 1;
    if plotCvAEquation
        clear CA; clear C; clear A; 
        Agrain = 0.001; 
        A = minA:Agrain:maxA;
        if lambdaC==0 && lambdaA==0 && gammaA==0.5
            top1 = rA^-betaC;
            top2inside = .5./(1-A);
            top2 = reallog(top2inside).^rB; 
            CA = 1-exp(-(top1.*top2));
        else
            top1 = rA^-betaC;
            top2inside = (1-gammaA-lambdaA)./(1-A-lambdaA); 
            top2 = (log(top2inside)).^rB;
            CA = gammaC + (1-gammaC-lambdaC) * (1-exp(-(top1.*top2))); 
        end
        plot(A,CA,':','Color',p.style.attColors(iAtt,:),'LineWidth',2)
        
        % Draw patch curved
        aPlot = area(A,CA);
        aPlot.FaceColor = p.style.labelColors; 
        % aPlot.FaceAlpha = 0.7; % alpha 
        aPlot.EdgeColor = p.style.labelColors; % [0 0 0]; 
        uistack(aPlot,'bottom')

        % Calculate AUC 
        Q(iAtt,:) = cumtrapz(A,CA); 
        AUC(iAtt) = real(Q(iAtt,end)); 
    end 
end

uistack(sPlot,"top")

%% Save analysis 
% a.rangeC = rangeC; 
% a.rangeA = rangeA;
% a.minC = minC; 
% a.maxC = maxC; 
% a.minA = minA; 
% a.maxA = maxA; 
aGroup.Cfit = Cfit; 
aGroup.Afit = Afit; 
aGroup.CA = CA; 
aGroup.AUC = AUC; 

