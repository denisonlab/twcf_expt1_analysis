function a = twcf_plotFitGroupAbsent(var,varName,varShortName,a,p,taskType,s)
% a = twcf_plotFit(var,varName,a,p,taskType,figureType)
% Plots var (e.g. A(x), C(x) by stim strength (x) by attention condition 
% And fits psychometric function
% 
% Inputs
%   a analysis structure
%   must contain task type 
%   'detection' or 'discrimination'
% Output 
% Expt 1: cue tex det 

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
        stimLevels = log(a.uniqueLineLengths_inDeg_postcue); % log space (log10) 
    case 'weibull'
        stimLevels = a.uniqueLineLengths_inDeg_postcue; % linear space <-- here subtract 1 pix? 
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
% subplot (1,1,1)
hold on

switch taskType
    case 'detection'
        ylims = [0 1]; 
    case 'discrimination'
        ylims = [0.4 1]; 
        yline(0.5,'--')
end
ylim(ylims)
yticks(0:0.25:1)
figureStyle
ylabel(sprintf('p(%s)', varShortName),'FontSize',p.style.textAxisSize)

switch p.fit.PFtype 
    case 'gumbel'
        xlabel('log_{10} (line length) at cued loc (°)','FontSize',p.style.textAxisSize)
        xlim([stimLevels(1)-p.style.xBuffer stimLevels(end)+p.style.xBuffer])
    case 'weibull'
        set(gca,'Xscale','log')
        xlabel('Line length (°)','FontSize',p.style.textAxisSize)
        xlim([stimLevels(1)-p.style.xBuffer/10 stimLevels(end)+p.style.xBuffer/10])
end
xl = get(gca, 'XLim');

xticks(stimLevels([1,2,4,6,7]))
xtickformat('%.2f')
xtickangle(p.style.xtickangle)

xAX = get(gca,'XAxis');  
yAX = get(gca,'YAxis');  
xlab = get(gca,'XLabel');
ylab = get(gca,'YLabel');

set(xAX,'FontSize', p.style.textTickSize)
set(yAX,'FontSize', p.style.textTickSize)
set(xlab,'FontSize', p.style.textAxisSize)
set(ylab,'FontSize', p.style.textAxisSize)

NumPos = []; OutOfNum = [];
OutOfNum = a.outOfNum'; % sum(a.stimPresent{i,iAtt},'omitnan'); % Number of trials at each entry of 'StimLevels'
for iAtt = 1:numel(a.attHeaders)
    for iL = 1:numel(stimLevels)
        x = stimLevels(iL);
        y = var(iL,iAtt);
        a.(varName).yMean(iL,iAtt) = y; 
        % Format values for Palamedes PF fit
        NumPos(iAtt,iL) = round(y*OutOfNum(iAtt,iL)); % sum(var{i,iAtt},'omitnan');
    end
    switch taskType
        case 'detection',      searchGrid.gamma = 0; % guess-rate
        case 'discrimination', searchGrid.gamma = 0.5;
    end
    % Fit a psychometric function to data using a Maximum Likelihood criterion
    [paramsFitted, LL, exitflag] = PAL_PFML_Fit(stimLevels, NumPos(iAtt,:), OutOfNum(iAtt,:), searchGrid, paramsFree, PF);
    fit = PF(paramsFitted, stimLevelsFine);

    l(iAtt) = plot(stimLevelsFine,fit,'Color',p.style.attColors(iAtt,:),'lineStyle','--','linewidth',p.style.fitLineWidth/2);
    
    % fitted alpha (threshold) 
    % xline(paramsFitted(1),'--','Color',p.style.attColors(iAtt,:),'linewidth',1.5) % threshold
   
%     findThreshold = paramsFitted(1) < stimLevelsFine;
%     thresholdIdx = find(findThreshold, 1, 'first');
% 
%     x = [paramsFitted(1) paramsFitted(1)];   % adjust length and location of arrow 
%     y = [fit(thresholdIdx) ylims(1)];
% 
%     [y_norm] = y_to_norm_v2(y(1),y(2));
%     [x_norm] = x_to_norm_v2(x(1),x(2));
% 
%     annotation('textarrow',x_norm,y_norm,'FontSize',13,'Linewidth',2,'Color',p.style.attColors(iAtt,:),'LineStyle',':')

    for iL = 1:numel(stimLevels)
        x = stimLevels(iL);
        y = var(iL,iAtt);
        scatter(x,y,p.style.sz*0.8,'x','MarkerFaceColor','w','MarkerFaceAlpha',1,'MarkerEdgeColor',p.style.attColors(iAtt,:),'LineWidth',2);
    end
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

% paramText = sprintf('\x03b1 = %0.2f, %0.2f, %0.2f; \x03b2 = %0.2f, %0.2f, %0.2f\n\x03b3 = %0.2f, %0.2f, %0.2f; \x03bb = %0.2f, %0.2f, %0.2f',...
%     a.(varName).paramsFitted(1,1),a.(varName).paramsFitted(2,1),a.(varName).paramsFitted(3,1),...
%     a.(varName).paramsFitted(1,2),a.(varName).paramsFitted(2,2),a.(varName).paramsFitted(3,2),...
%     a.(varName).paramsFitted(1,3),a.(varName).paramsFitted(2,3),a.(varName).paramsFitted(3,3),...
%     a.(varName).paramsFitted(1,4),a.(varName).paramsFitted(2,4),a.(varName).paramsFitted(3,4)); 
% title(paramText) 

%% Save analysis 
a.(varName).var = var; 
a.(varName).varName = varName; 
a.(varName).varShortName = varShortName; 
a.(varName).figureType = s.figureType; 
a.(varName).taskType = taskType; 

a.(varName).stimLevels = stimLevels; 
a.(varName).stimLevelsFine = stimLevelsFine; 
a.(varName).attHeaders = a.attHeaders; 

a.(varName).PF = PF; 
a.(varName).searchGrid = searchGrid; 
a.(varName).paramValueNames = {'threshold','slope','guess-rate','lapse-rate'};
a.(varName).paramsFree = paramsFree; 

a.(varName).xl = xl; 

