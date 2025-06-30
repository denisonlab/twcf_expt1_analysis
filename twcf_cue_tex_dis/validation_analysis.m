function validation_analysis(dataFile, baseDir)

% Validation is all neutral cue 

%% Check inputs
if nargin<1
    error('Must provide data file');
end

%% load data
saveFigs       = 1; 
exptName       = 'cued texture discrimination';
exptShortName  = 'twcf_cue_tex_dis'; 
site           = 'BU';
dataFolder     = 'validation';
taskType       = 'discrimination'; 

% * * * path * * *
% currentFolder = '/Users/karen/Documents/GitHub/twcf_expt1_data_BU'; % lab desk
% currentFolder = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_data_BU'; % laptop
currentFolder = sprintf('%s/twcf_expt1_data_BU',baseDir); 
% * * * * * *

dataFileSplitName = strsplit(dataFile,'_');
subjectID  = dataFileSplitName{6}; % 'S0006';
date       = dataFileSplitName{8}; % '20220614_153756'; 
time       = dataFileSplitName{9}; % '20220614_153756'; \

subjectFolder = sprintf('%s/%s/%s/%s',currentFolder,exptShortName,subjectID,dataFolder);  
load(sprintf('%s/%s.mat',subjectFolder,dataFile))

addpath (genpath(pwd))

if saveFigs 
    figFolder = sprintf('%s/figs',subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% Calculate accuracy by line length 
a.uniqueLineLengths_inDeg_postcue = unique(data.lineLength_inDeg_postcue); % 7 lengths 
% no accuracy measure for orientation discrim when stim is absent = -1 
% no keyboard input = -2

% break down discrim and detect measures by line length and att conds
a.attConds = unique(data.cueValidity); % 1 --> valid, 0 --> neutral, -1 --> invalid
a.attHeaders = {'neutral'};
for i = 1:numel(a.uniqueLineLengths_inDeg_postcue)
    % all
    a.correctDis{i} = data.correctDis(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.correctDis~=-2);
    a.ratingDis{i} = data.ratingDis(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1);
    a.cueValidity{i} = data.cueValidity(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1);
    
    for iAtt = 1:numel(a.attHeaders)
        a.correctDis_att{i,iAtt} = data.correctDis(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.correctDis~=-2 & data.cueValidity==a.attConds(iAtt));
        a.ratingDis_att{i,iAtt} = data.ratingDis(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.cueValidity==a.attConds(iAtt));
        a.cueValidity_att{i,iAtt} = data.cueValidity(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.cueValidity==a.attConds(iAtt));
        
        % Detection 
%         a.stimPresent{i,iAtt} = stimPresent(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%         a.sawFigure{i,iAtt} = sawFigure(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%         a.sawShape{i,iAtt} = sawShape(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%         
%         a.sawFigure_stimPresent{i,iAtt} = sawFigure_stimPresent(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%         a.sawShape_stimPresent{i,iAtt} = sawShape_stimPresent(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%     
%         a.sawFigure_stimAbsent{i,iAtt} = sawFigure_stimAbsent(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
%         a.sawShape_stimAbsent{i,iAtt} = sawShape_stimAbsent(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.goodTrial==1 & data.cueValidity==a.attConds(iAtt)); 
        
        % Discrimination 
        % 1-->vertical, 0-->horizontal
        a.stimOri{i,iAtt} = data.stimID_postcue(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.cueValidity==a.attConds(iAtt) & data.goodTrial==1);
        % 123-->H, 456-->V
        a.responseOri{i,iAtt} = data.respDis(data.lineLength_inDeg_postcue==a.uniqueLineLengths_inDeg_postcue(i) & data.cueValidity==a.attConds(iAtt) & data.goodTrial==1);
    end
end

%% Plot objective subjective 
% Figure styling
sz = 100;
alpha = 0.5; 
xBuffer = 0.1; 
attColors = {'k'}; % only neutral 

% Figure
figure
set(gcf,'Position',[100 100 700 400])
sgtitle(sprintf('expt: %s | site: %s \n subject: %s | validation', exptName, site, subjectID))
% p(correct) discrimination 
subplot 211
hold on
ylim([0 1])
xlim([log(a.uniqueLineLengths_inDeg_postcue(1))-xBuffer log(a.uniqueLineLengths_inDeg_postcue(end))+xBuffer])
yticks([0 0.5 1])
yline(0.5,'--')
variable = a.correctDis_att; 
NumPos = []; OutOfNum = [];
for iAtt = 1:numel(a.attHeaders)
    for i = 1:numel(a.uniqueLineLengths_inDeg_postcue)
        x = log(a.uniqueLineLengths_inDeg_postcue(i));
        y = variable{i,iAtt};
        y = mean(y,'omitnan'); 
        pCorrect(i,iAtt) = y; % save these values
        scatter(x,y,sz,'MarkerFaceColor',attColors{iAtt},'MarkerFaceAlpha',alpha,'MarkerEdgeColor','w');
        % Format values for Palamedes PF fit
        % NumPos(iAtt,i) = sum(variable{i,iAtt},'omitnan');
        % OutOfNum(iAtt,i) = sum(variable{i,iAtt},'omitnan'); % sum(a.stimPresent{i,iAtt},'omitnan'); % Number of trials at each entry of 'StimLevels'
    end
    switch taskType
        case 'detection',      searchGrid.gamma = 0; % guess-rate
        case 'discrimination', searchGrid.gamma = 0.5;
    end
%     [paramsValues LL exitflag] = PAL_PFML_Fit(stimLevels, NumPos(iAtt,:), OutOfNum(iAtt,:), searchGrid, paramsFree, PF);
%     Fit = PF(paramsValues, stimLevelsFine);
%     l(iAtt) = plot(stimLevelsFine,Fit,attColors{iAtt},'linewidth',1);
%     xline(paramsValues(1),'--','Color',attColors{iAtt}) % threshold
%     legend([l(1),l(2),l(3)],a.attHeaders,'location','southeast')
end
ylabel(sprintf('p(correct)\ndiscrimination'))
xlabel('log_{10} (line length) at cued loc (deg)')
xl = get(gca, 'XLim');
figureStyle
title(sprintf('avg # of trials per data point = %d (invalid), %d (neutral), %d (valid)',...
    0,...
    sum(data.goodTrial),...
    0))

% p(stronger), subjective rating 
subplot 212
hold on
ylim([0 1])
xlim([log(a.uniqueLineLengths_inDeg_postcue(1))-xBuffer log(a.uniqueLineLengths_inDeg_postcue(end))+xBuffer])
yticks([0 0.5 1])
yline(0.5,'--')
variable = a.ratingDis_att; 
NumPos = []; OutOfNum = [];
for iAtt = 1:numel(a.attHeaders)
    for i = 1:numel(a.uniqueLineLengths_inDeg_postcue)
        x = log(a.uniqueLineLengths_inDeg_postcue(i));
        y = variable{i,iAtt};
        y = mean(y,'omitnan'); 
        pStronger(i,iAtt) = y; % save these values
        scatter(x,y,sz,'MarkerFaceColor',attColors{iAtt},'MarkerFaceAlpha',alpha,'MarkerEdgeColor','w');
        % Format values for Palamedes PF fit
        % NumPos(iAtt,i) = sum(variable{i,iAtt},'omitnan');
        % OutOfNum(iAtt,i) = sum(variable{i,iAtt},'omitnan'); % sum(a.stimPresent{i,iAtt},'omitnan'); % Number of trials at each entry of 'StimLevels'
    end
    switch taskType
        case 'detection',      searchGrid.gamma = 0; % guess-rate
        case 'discrimination', searchGrid.gamma = 0.5;
    end
%     [paramsValues LL exitflag] = PAL_PFML_Fit(stimLevels, NumPos(iAtt,:), OutOfNum(iAtt,:), searchGrid, paramsFree, PF);
%     Fit = PF(paramsValues, stimLevelsFine);
%     l(iAtt) = plot(stimLevelsFine,Fit,attColors{iAtt},'linewidth',1);
%     xline(paramsValues(1),'--','Color',attColors{iAtt}) % threshold
%     legend([l(1),l(2),l(3)],a.attHeaders,'location','southeast')
end
ylabel(sprintf('p(stronger)\nsubjective'))
xlabel('log_{10} (line length) at cued loc (deg)')
xl = get(gca, 'XLim');
figureStyle
title(sprintf('avg # of trials per data point = %d (invalid), %d (neutral), %d (valid)',...
    0,...
    sum(data.goodTrial),...
    0))

if saveFigs
    figTitle = sprintf('%s_%s_%s_%s_%s_%s',exptShortName,site,subjectID,date,time,dataFolder);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

