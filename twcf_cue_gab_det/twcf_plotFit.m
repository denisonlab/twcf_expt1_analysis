function a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s)
% a = twcf_plotFit(var,varName,a,p,taskType,figureType)
% Plots var (e.g. A(x), C(x) by stim strength (x) by attention condition 
% And fits psychometric function
% 1.3 cue gab det
% 
% Inputs
%   a analysis structure
%   must contain task type 
%   'detection' or 'discrimination'
% Output 

%% Check inputs
if nargin < 1 
    error('Must provide variable to plot')
end
if nargin < 2
    varName = 'variable name';
end
if nargin < 3 
    varShortName = 'var';
end

%% Palamedes fit parameters setup (move to params file) 
PF = p.fit.PF; % @PAL_Gumbel; % psychometric function type to fit, @PAL_Logistic, @PAL_Gumbel (ie log-Weibull), @PAL_Weibull

% Stimulus intensities and search grid 
switch p.fit.PFtype 
    case 'gumbel'
        a.uniqueContrasts_postcue(1) = a.uniqueContrast_postcue+0.01; % log of 0 is Inf, so set a correction to be 0.001
        stimLevels = log(a.uniqueContrasts_postcue); 
    case {'weibull','nakarushton'}
        stimLevels = a.uniqueContrasts_postcue; 
end
stimLevelsFine = min(stimLevels):(max(stimLevels)-min(stimLevels))/1000:max(stimLevels);

if contains(varName,'A') % accuracy 
    paramsFree        = p.fit.paramsFreeA; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed 
elseif contains(varName,'C') % consciousness 
    paramsFree        = p.fit.paramsFreeC;
end
searchGrid.alpha  = p.fit.searchGrid.alpha; % threshold 
searchGrid.beta   = p.fit.searchGrid.beta; % slope (adjust upper bound) 
searchGrid.lambda = p.fit.searchGrid.lambda; % lapse-rate 

%% Plot
subplot (4,1,p.iPlot)
hold on
switch taskType
    case 'detection'
        ylim([0 1])
    case 'discrimination'
        ylim([0.3 1])
end
yticks([0 0.5 1])
yline(0.5,'--')
figureStyle
ylabel(sprintf('%s\np(%s | \n %s)', varName, varShortName, s.gratingType))

NumPos = []; OutOfNum = [];
for iAtt = 1:numel(a.attHeaders)
    for i = 1:numel(stimLevels)
        x = stimLevels(i);
        y = var{i,iAtt};
        y = y(~isnan(y)); 
        y = mean(y,'omitnan');
        scatter(x,y,p.style.sz,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha,'MarkerEdgeColor','w');
        a.(varName).yMean(i,iAtt) = y; 
        % Format values for Palamedes PF fit
        NumPos(iAtt,i) = sum(var{i,iAtt},'omitnan');
        if i==1 
            OutOfNum(iAtt,i) = sum(~a.stimPresent{i,iAtt},'omitnan'); 
        else
        OutOfNum(iAtt,i) = sum(a.stimPresent{i,iAtt},'omitnan'); % Number of trials at each entry of 'StimLevels', change to only good trials?    
        end
    end
    switch taskType
        case 'detection',      searchGrid.gamma = 0:.1:1; % guess-rate
        case 'discrimination', searchGrid.gamma = 0.5;
    end
    % Fit a psychometric function to data using a Maximum Likelihood criterion
    [paramsFitted, LL, exitflag] = PAL_PFML_Fit(stimLevels, NumPos(iAtt,:), OutOfNum(iAtt,:), searchGrid, paramsFree, PF);
    fit = PF(paramsFitted, stimLevelsFine);
    l(iAtt) = plot(stimLevelsFine,fit,'Color',p.style.attColors(iAtt,:),'linewidth',1);
    xline(paramsFitted(1),'--','Color',p.style.attColors(iAtt,:),'linewidth',1.5) % threshold

    a.(varName).paramsFitted(iAtt,:) = paramsFitted; 
    a.(varName).LL(iAtt,:) = LL; 
    a.(varName).exitflag(iAtt,:) = exitflag; 
    a.(varName).fit(iAtt,:) = fit; 
end
a.(varName).NumPos = NumPos;
a.(varName).OutOfNum = OutOfNum;

% if p.iPlot==1
%     legend([l(1),l(2),l(3)],a.attHeaders,'location','southeast')
%     title(sprintf('avg # of trials per data point = %d (invalid), %d (neutral), %d (valid)',...
%     round(mean(OutOfNum(find(ismember(a.attHeaders, 'invalid')),:))),...
%     round(mean(OutOfNum(find(ismember(a.attHeaders, 'neutral')),:))),...
%     round(mean(OutOfNum(find(ismember(a.attHeaders, 'valid')),:)))))
% end

paramText = sprintf('\x03b1 = %0.2f, %0.2f, %0.2f; \x03b2 = %0.2f, %0.2f, %0.2f\n\x03b3 = %0.2f, %0.2f, %0.2f; \x03bb = %0.2f, %0.2f, %0.2f',...
    a.(varName).paramsFitted(1,1),a.(varName).paramsFitted(2,1),a.(varName).paramsFitted(3,1),...
    a.(varName).paramsFitted(1,2),a.(varName).paramsFitted(2,2),a.(varName).paramsFitted(3,2),...
    a.(varName).paramsFitted(1,3),a.(varName).paramsFitted(2,3),a.(varName).paramsFitted(3,3),...
    a.(varName).paramsFitted(1,4),a.(varName).paramsFitted(2,4),a.(varName).paramsFitted(3,4)); 

if p.iPlot==1
    title(sprintf('Fitted parameters (invalid, neutral, valid)\n%s',paramText)); 
else
    title(paramText) 
end

clear xticklabel
for iStim = 1:numel(stimLevels)
    xticklabel{iStim,:} = sprintf('%0.3f',stimLevels(iStim)); 
end
xticks(stimLevels)
xticklabels(xticklabel)

switch p.fit.PFtype 
    case 'gumbel'
        xlabel('log_{10} (line length) at cued loc (deg)','FontSize',p.style.textAxisSize)
        xlim([stimLevels(1)-p.style.xBuffer stimLevels(end)+p.style.xBuffer])
    case 'weibull'
        switch s.exptShortName
            case {'twcf_cue_tex_det','twcf_cue_tex_dis'}
                set(gca,'Xscale','log')
                xlabel('Line length (deg)','FontSize',p.style.textAxisSize)
                xlim([stimLevels(1)-p.style.xBuffer/10 stimLevels(end)+p.style.xBuffer/10])
            case {'twcf_cue_gab_det','twcf_cue_gab_dis'}
                xlabel('Contrast','FontSize',p.style.textAxisSize)
                xlim([stimLevels(1)-p.style.xBufferContrast stimLevels(end)+p.style.xBufferContrast])
        end

end

xl = get(gca, 'XLim');

%% Save analysis 
a.(varName).var = var; 
a.(varName).varName = varName; 
a.(varName).varShortName = varShortName; 
a.(varName).gratingType = s.gratingType; 
a.(varName).taskType = taskType; 

a.(varName).stimLevels = stimLevels; 
a.(varName).stimLevelsFine = stimLevelsFine; 
a.(varName).attHeaders = a.attHeaders; 

a.(varName).PF = PF; 
a.(varName).searchGrid = searchGrid; 
a.(varName).paramValueNames = {'threshold','slope','guess-rate','lapse-rate'};
a.(varName).paramsFree = paramsFree; 

a.(varName).xl = xl; 







