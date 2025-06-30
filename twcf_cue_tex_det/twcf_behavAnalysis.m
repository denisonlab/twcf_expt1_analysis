function [a,PSI] = twcf_behavAnalysis(dataAll,subjectIDs,subjectID,s)
% TWCF 1.1 cue tex det 
% [a] = behav_analysis(dataAll,subjectIDs,subjectID)

% runBehav
% Inputs 
%   dataAll is structure of concatenated session data 
%   subjectIDs is cell array of subject IDs 
%   subjectID is string of subject of interest 
%   s is setup structure
% Outputs
%   a is analysis structure 

%% Do analysis on figure present or figure absent trials (move to setup) 
s.figureType     = 'figure present';  % figure absent, figure present 

%% Setup 
addpath(genpath(pwd))

% Get figure and psychometric function fitting parameters 
p = twcf_analysisParams; 

%% Load data 
% % * * * change here * * *
% subjectID = 'S0004'; % subject of interest
% % * * * * * * * * * * * *
subjectIdx = find(ismember(subjectIDs, subjectID));
clear data 
data = dataAll(subjectIdx);
s.subjectFolder = sprintf('%s/%s/%s/%s',s.dataFolder,s.exptShortName,subjectID,s.exptFolder);

% Make figure directory 
if s.saveFigs 
    figFolder = sprintf('%s/figs',s.subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% Set different conditions of interest and condition names 
% C
cVarNames      = {'sawFigure','sawShape','sawFigureNS'}; 
cFieldNames    = {'C1','C2','CX'}; 
cVarLabelNames = {'p(saw figure)','p(saw shape)','p(saw figure not shape)'}; 

% A 
aVarNames      = {'correctDis'};
aFieldNames    = {'A'}; 
aVarLabelNames = {'p(correct discrimination)'}; 

% All var names 
varNames       = {'A','C1','C2','CX'};

% Attention 
a.attConds   = unique(data.cueValidity); % -1 --> invalid, 1 --> valid, 0 --> neutral
a.attHeaders = {'invalid','neutral','valid'};

%% Get indices stim present, saw figure, saw shape, saw figure | saw shape from response buttons 
idx.stimPresent = logical(data.stimID_postcue(:)>0)';
idx.stimAbsent  = logical(data.stimID_postcue(:)==0)';

idx.sawFigure           = double(logical(data.resp==1 | data.resp==2 | data.resp==5 | data.resp==6)); % if saw shape, also saw figure 
idx.sawFigureNS         = double(logical(data.resp==2 | data.resp==5)); % saw figure and not shape 
idx.sawShape            = double(logical(data.resp==1 | data.resp==6)); 

switch s.figureType 
    case 'figure present'
        idxCase = ~idx.stimPresent; 
    case 'figure absent'
        idxCase = ~idx.stimAbsent; 
end
idx.sawFigure_present = idx.sawFigure; 
idx.sawFigure_absent  = idx.sawFigure; 

idx.sawFigureNS_present = idx.sawFigureNS; 
idx.sawFigureNS_absent  = idx.sawFigureNS; 

idx.sawShape_present = idx.sawShape; 
idx.sawShape_absent  = idx.sawShape; 

idx.sawFigure_present(~idx.stimPresent)         = NaN;
idx.sawFigureNS_present(~idx.stimPresent)       = NaN;
idx.sawShape_present(~idx.stimPresent)          = NaN;

idx.sawFigure_absent(~idx.stimAbsent)           = NaN;
idx.sawFigureNS_absent(~idx.stimAbsent)         = NaN;
idx.sawShape_absent(~idx.stimAbsent)            = NaN;

%% Break down data by line length and attention conditions 
a.uniqueLineLengths_inDeg_postcue = unique(data.lineLength_inDeg_postcue); % 7 stimulus strengths 
for i = 1:numel(a.uniqueLineLengths_inDeg_postcue)
    if a.uniqueLineLengths_inDeg_postcue(i) > 3
        a.uniqueLineLengths_inDeg_postcue(i) = 3; % correction for if threshold surpasses maximum possible line length
    end
end

a.log10_uniqueLineLengths_inDeg_postcue = log10(unique(data.lineLength_inDeg_postcue));
% no accuracy measure for orientation discrim when stim is absent = -1; no keyboard input = -2

idx.good = data.goodTrial==1;
for i = 1:numel(a.uniqueLineLengths_inDeg_postcue)
    % all trials 
    lines =         data.lineLength_inDeg_postcue; 
    lines(lines>3) = 3; 
    idx.line = lines==a.uniqueLineLengths_inDeg_postcue(i);

    a.correctDis_allTrials{i}     = data.correctDis(idx.line & idx.good & data.correctDis~=-1 & idx.stimPresent);
    a.correctDet_allTrials{i}     = data.correctDet(idx.line & idx.good & idx.stimPresent);
    a.cueValidity_allTrials{i}    = data.cueValidity(idx.line & idx.good & idx.stimPresent);
    
    a.sawFigure_allTrials{i}      = idx.sawFigure(idx.line & idx.good & idx.stimPresent); 
    a.sawFigureNS_allTrials{i}    = idx.sawFigureNS(idx.line & idx.good & idx.stimPresent); 
    a.sawShape_allTrials{i}       = idx.sawShape(idx.line & idx.good & idx.stimPresent); 
    
    % trials by attention 
    for iAtt = 1:numel(a.attHeaders)
        idx.att = data.cueValidity==a.attConds(iAtt);
        a.stimPresent{i,iAtt}     = idx.stimPresent(idx.line & idx.good & idx.att & idx.stimPresent); 
        a.goodTrial{i,iAtt}       = data.goodTrial(idx.line & idx.att & idx.stimPresent);
        a.cueValidity{i,iAtt}     = data.cueValidity(idx.line & idx.good & idx.att & idx.stimPresent);

        a.correctDis{i,iAtt}      = data.correctDis(idx.line & idx.good & data.correctDis~=-1 & idx.att & idx.stimPresent);
        a.correctDet{i,iAtt}      = data.correctDet(idx.line & idx.good & idx.att & idx.stimPresent);
          
        % Detection (figure present) 
        a.sawFigure{i,iAtt}       = idx.sawFigure_present(idx.line & idx.good & idx.att & idx.stimPresent); 
        a.sawFigureNS{i,iAtt}     = idx.sawFigureNS_present(idx.line & idx.good & idx.att & idx.stimPresent); 
        a.sawShape{i,iAtt}        = idx.sawShape_present(idx.line & idx.good & idx.att & idx.stimPresent); 
        
        % Detection (all trials figure present and absent) to calculate d'
        % and c measures 
        a.sawFigure_absent{i,iAtt}       = idx.sawFigure_absent(idx.line & idx.good & idx.att & idx.stimAbsent); 
        a.sawFigureNS_absent{i,iAtt}     = idx.sawFigureNS_absent(idx.line & idx.good & idx.att & idx.stimAbsent); 
        a.sawShape_absent{i,iAtt}        = idx.sawShape_absent(idx.line & idx.good & idx.att & idx.stimAbsent); 
        
        % Detection (both figure present and absent trials) 
        % to calculate d' and c measures 
        a.sawFigure_presentAbsent{i,iAtt}       = idx.sawFigure(idx.line & idx.good & idx.att); 
        a.sawFigureNS_presentAbsent{i,iAtt}     = idx.sawFigureNS(idx.line & idx.good & idx.att); 
        a.sawShape_presentAbsent{i,iAtt}        = idx.sawShape(idx.line & idx.good & idx.att);

        % Discrimination 
        a.stimOri{i,iAtt}        = data.stimID_postcue(idx.line & idx.att & idx.good & idx.stimPresent); % 0-->absent, 1-->vertical, 2-->horizontal
        a.responseOri{i,iAtt}    = data.respDis(idx.line& idx.att &idx.good & idx.stimPresent); % 123-->H, 456-->V
        
        % For calcluating discrimiantion SDT measures 
        a.dis_correct_V{i,iAtt} = data.correctDis(idx.line& idx.att &idx.good & idx.stimPresent & data.stimID_postcue==1); 
        a.dis_correct_H{i,iAtt} = data.correctDis(idx.line& idx.att &idx.good & idx.stimPresent & data.stimID_postcue==2); 

        % RT 
        a.RT{i,iAtt} = data.RT(idx.line & idx.good & idx.att & idx.stimPresent); 
    end
end

for i = 1:numel(cFieldNames)
    n.(cVarNames{i}) = sum(cellfun(@sum, a.(cVarNames{i})),'all','omitnan'); % get detection resp type counts 
end

%% Calculate SDT measures
clear sdt

% Detection SDT --> recalculate collapsed across line lengths for absent 
for iC = 1:numel(cFieldNames)
    for iL = 1:numel(a.uniqueLineLengths_inDeg_postcue)
        for iAtt = 1:numel(a.attHeaders)
            condFA = []; nfa = []; nnoise = [];
            condH = []; nh = []; nsignal = [];
            
            condName = sprintf('%s_absent',cVarNames{iC}); 
            condFA = a.(condName){iL,iAtt};
            nfa = numel(find(condFA==1));
            nnoise = numel(condFA); 
            
            condName = sprintf('%s',cVarNames{iC}); 
            condH = a.(condName){iL,iAtt};
            nh = numel(find(condH==1));
            nsignal = numel(condH);
            
            loglinearCorrection = 1; 
            [sdt.(condName).dprimeDetect(iL,iAtt), sdt.(condName).criterionDetect(iL,iAtt)] = kt_dprime(nh,nfa,nsignal,nnoise,loglinearCorrection);
            
            sdt.(condName).det_nh(iL,iAtt) = nh;
            sdt.(condName).det_nfa(iL,iAtt) = nfa;
            sdt.(condName).det_nsignal(iL,iAtt) = nsignal;
            sdt.(condName).det_nnoise(iL,iAtt) = nnoise;
            sdt.(condName).det_ph(iL,iAtt) = nh/nsignal;
            sdt.(condName).det_pfa(iL,iAtt) = nfa/nnoise;
        end
    end
end

% Discrimination SDT 
for iL = 1:numel(a.uniqueLineLengths_inDeg_postcue)
    for iAtt = 1:numel(a.attHeaders)
        condFA = []; nfa = []; nnoise = [];
        condH = []; nh = []; nsignal = [];
  
        condName = 'dis_correct_H'; % Take one orientation to be signal, other as noise 
        condFA = a.(condName){iL,iAtt};
        nfa = numel(find(condFA==0));
        nnoise = numel(condFA);
        
        condName = 'dis_correct_V'; % Take one orientation to be signal, other as noise 
        condH = a.(condName){iL,iAtt};
        nh = numel(find(condH==1));
        nsignal = numel(condH);
        
        loglinearcorrection = 1; 
        [sdt.(condName).dprimeDis(iL,iAtt), sdt.(condName).criterionDis(iL,iAtt)] = kt_dprime(nh,nfa,nsignal,nnoise,loglinearcorrection);

        sdt.(condName).dis_nh(iL,iAtt) = nh;
        sdt.(condName).dis_nfa(iL,iAtt) = nfa;
        sdt.(condName).dis_nsignal(iL,iAtt) = nsignal;
        sdt.(condName).dis_nnoise(iL,iAtt) = nnoise;
        sdt.(condName).dis_ph(iL,iAtt) = nh/nsignal;
        sdt.(condName).dis_pfa(iL,iAtt) = nfa/nnoise;
    end
end

% Save 
a.sdt = sdt; 

%% Setup figure 
figure
set(gcf,'Position',[0 0 500 p.style.figHeight])
sgtitle(sprintf('expt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 0; 

%% Discrimination, A, p(correct)
p.iPlot = p.iPlot+1; 
var = a.correctDis; 
varName = 'A'; % discriminationpCorrect
varShortName = 'correct';
taskType = 'discrimination';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Detection, C1, p(saw figure)
p.iPlot = p.iPlot + 1; 
var = a.sawFigure; 
varName = 'C1'; % 
varShortName = 'saw figure';
taskType = 'detection';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Detection, C2, p(saw shape)
p.iPlot = p.iPlot+1; 
var = a.sawShape; 
varName = 'C2'; 
varShortName = 'saw shape';
taskType = 'detection';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

% annotation = sprintf('%d of %d trials (%0.2f)',n.(cVarNames{2}),n.(cVarNames{1}),n.(cVarNames{2})/n.(cVarNames{1}));
% t = title(annotation);
% t.FontWeight = 'normal'; 

%% Detection, CX, p(saw figure and not saw shape)
p.iPlot = p.iPlot+1; 
var = a.sawFigureNS; 
varName = 'CX'; 
varShortName = 'saw figure and not shape';
taskType = 'detection';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

% annotation = sprintf('%d of %d trials (%0.2f)',n.(cVarNames{3}),n.(cVarNames{1}),n.(cVarNames{3})/n.(cVarNames{1}));
% t = title(annotation);
% t.FontWeight = 'normal'; 

%% Save figure with detection discrimination subplots 
if s.saveFigs 
    figTitle = sprintf('%s_%s_%s_%s_CAbyXAtt',s.exptFolder,s.site,subjectID,datestr(now,'yymmdd')); 
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType)) 
end

%% Plots fitted params (subject-level)
figure
set(gcf,'Position',[0 0 400 p.style.figHeight])
sgtitle(sprintf('Fitted params\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 1;

for i = 1:numel(varNames)
    for iP = 1:numel(a.(varNames{i}).paramValueNames)
        subplot (numel(varNames),numel(a.(varNames{i}).paramValueNames),p.iPlot)
        hold on
        figureStyle
        clear x; clear y;
        y = a.(varNames{i}).paramsFitted;
        xlim([iP-1 iP+1])
        xticks(iP)
        xticklabels(a.(varNames{i}).paramValueNames{iP})
        ylabel(sprintf(varNames{i}))   
        for iAtt = 1:numel(a.(varNames{i}).attHeaders)
            scatter(iP,y(iAtt,iP),p.style.sz,'filled','MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha,'MarkerEdgeColor','w')
        end
        p.iPlot = p.iPlot+1;
    end

end

if s.saveFigs
    figTitle = sprintf('%s_%s_%s_%s_fittedParams',s.exptFolder,s.site,subjectID,datestr(now,'yymmdd'));
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType))
end

%% Plots CvA
% Testing performance-relative subjective inflation (PSI) 
figure
set(gcf,'Position',[0 0 350 p.style.figHeight])
sgtitle(sprintf('Testing PSI - Metaperceptual function (CvA)\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 1; 

varA = a.A; 
for iC = 1:numel(cFieldNames) 
    varC = a.(cFieldNames{iC}); 
    PSI.(cFieldNames{iC}) = twcf_plotCvA(varA,varC,p,a);
    p.iPlot = p.iPlot+1;
end

if s.saveFigs
    figTitle = sprintf('%s_%s_%s_%s_CvA',s.exptFolder,s.site,subjectID,datestr(now,'yymmdd'));
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType))
end

%% Plots AUC 
figure
set(gcf,'Position',[0 0 350 p.style.figHeight])
sgtitle(sprintf('Testing PSI - AUC of CvA\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 1; 

% Plots CvA
varA = a.A;
for iC = 1:numel(cFieldNames)
    varC = a.(cFieldNames{iC});
    % Plot AUC of CvA
    subplot (3,1,p.iPlot)
    hold on
    figureStyle
    axis square
    for iAtt = 1:numel(a.attHeaders)
        clear x; clear y;
        x = [1 2 3];
        y = PSI.(cFieldNames{iC}).AUC;
        scatter(x(iAtt),y(iAtt),p.style.sz,'filled','MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha,'MarkerEdgeColor','w')
    end
    xlim([0 4])
    ylabel(sprintf('AUC of %svA',cFieldNames{iC}))
    xticks([1,2,3])
    xticklabels({'I','N','V'})
    p.iPlot = p.iPlot+1;
end

if s.saveFigs
    figTitle = sprintf('%s_%s_%s_%s_AUC',s.exptFolder,s.site,subjectID,datestr(now,'yymmdd'));
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType ))
end

%% Save analyses 
% General data parameters 
a.exptName         = s.exptName; 
a.exptShortName    = s.exptShortName;
a.site             = s.site; 
a.exptFolder       = s.exptFolder; 
a.figureType       = s.figureType; 
a.dataFolder       = s.dataFolder; 
a.subjectFolder    = s.subjectFolder; 
a.subjectIdx       = subjectIdx; 
a.subjectData      = data; 



