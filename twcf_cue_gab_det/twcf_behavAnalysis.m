function [a,PSI] = twcf_behavAnalysis(dataAll, subjectIDs, subjectID, s, plotPSI)
% [a] = behav_analysis(dataAll,subjectIDs,subjectID)
% 1.3 cue gab det 

% runBehav
% Inputs 
%   dataAll is structure of concatenated session data 
%   subjectIDs is cell array of subject IDs 
%   subjectID is string of subject of interest 
%   s is setup structure
% Outputs
%   a is analysis structure 

%% Do analysis on figure present or figure absent trials (move to setup) 
s.gratingType     = 'grating present';  % 'grating present', 'grating absent'

%% Setup 
addpath(genpath(pwd))

% Get figure and psychometric function fitting parameters 
p = twcf_analysisParams; 

%% Load data 
subjectIdx = find(ismember(subjectIDs, subjectID));
clear data 
data = dataAll(subjectIdx);
s.subjectFolder = sprintf('%s/%s/corrected/%s/%s',s.dataFolder,s.exptShortName,subjectID,s.exptFolder);

% Make figure directory 
if s.saveFigs % change when using UCI+BU data to unlocked folder 
    figFolder = sprintf('%s/figs',s.subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% Set different conditions of interest and condition names 
% C
cVarNames      = {'sawGrating','sawOri','sawGratingNO'}; 
cFieldNames    = {'C1','C2','CX'}; 
cVarLabelNames = {'p(saw grating)','p(saw ori)','p(saw grating not ori)'}; 

% A 
aVarNames      = {'correctDis'};
aFieldNames    = {'A'}; 
aVarLabelNames = {'p(correct discrimination)'}; 

% All var names 
varNames       = {'A','C1','C2','CX'};

% Attention 
a.attConds   = unique(data.cueValidity); % -1 --> invalid, 1 --> valid, 0 --> neutral;
a.attHeaders = {'invalid','neutral','valid'};

%% Get indices stim present, saw figure, saw shape, saw figure | saw shape from response buttons 
idx.stimPresent = logical(data.stimID_postcue(:)>0)';
idx.stimAbsent  = logical(data.stimID_postcue(:)==0)';

idx.sawGrating           = double(logical(data.resp==1 | data.resp==2 | data.resp==5 | data.resp==6)); % if saw shape, also saw figure 
idx.sawGratingNO         = double(logical(data.resp==2 | data.resp==5)); % saw figure and not shape 
idx.sawOri               = double(logical(data.resp==1 | data.resp==6)); 

switch s.gratingType 
    case 'grating present'
        idxCase = ~idx.stimPresent; 
    case 'grating absent'
        idxCase = idx.stimPresent; 
end
idx.sawGrating_present = idx.sawGrating; 
idx.sawGrating_absent  = idx.sawGrating; 

idx.sawGratingNO_present = idx.sawGratingNO; 
idx.sawGratingNO_absent  = idx.sawGratingNO; 

idx.sawOri_present = idx.sawOri; 
idx.sawOri_absent  = idx.sawOri; 

idx.sawGrating_present(~idx.stimPresent)         = NaN;
idx.sawGratingNO_present(~idx.stimPresent)       = NaN;
idx.sawOri_present(~idx.stimPresent)             = NaN;

idx.sawGrating_absent(~idx.stimAbsent)           = NaN;
idx.sawGratingNO_absent(~idx.stimAbsent)         = NaN;
idx.sawOri_absent(~idx.stimAbsent)               = NaN;

%% Calculate precue contrast
for iTrial = 1:size(data.contrastID_all,2)
    if data.precueLoc(iTrial)<1
        data.contrast_precue(iTrial) = NaN; 
        data.stimID_precue(iTrial) = NaN; 
    else
        data.contrast_precue(iTrial) = data.contrastID_all(data.precueLoc(iTrial),iTrial);
        data.stimID_precue(iTrial) = data.stimID_all(data.precueLoc(iTrial),iTrial); 
    end

    if isnan(data.respDet(iTrial)) || isnan(data.contrast_precue(iTrial))
        data.correctDet_conditionedOnPrecue(iTrial) = NaN;
    elseif (data.respDet(iTrial)==1 && data.contrast_precue(iTrial)>0) || data.respDet(iTrial)==0 && data.contrast_precue(iTrial)==0
        data.correctDet_conditionedOnPrecue(iTrial) = 1; 
    else
        data.correctDet_conditionedOnPrecue(iTrial) = 0; 
    end
end

% data.correctDet_conditionedOnPrecue = double (data.respDet && data.contrast_precue>0);
% data.correctDet_conditionedOnPrecue(data.contrast_precue<=0) = NaN; 

data.correctDis_conditionedOnPrecue = double(data.respDis==data.stimID_precue); 
data.correctDis_conditionedOnPrecue(data.stimID_precue==0) = NaN; 

%% Calculate distractor contrasts 
for iTrial = 1:size(data.contrastID_all,2)
    data.contrast_distractors(:,iTrial) = data.contrastID_all(:,iTrial);
    data.contrast_distractors(data.postcueLoc(iTrial),iTrial) = NaN; 
end

%% Break down data by contrast and attention conditions 
a.uniqueContrasts_postcue = unique(data.contrast_postcue); % 7 stimulus strengths 
a.log10_uniqueContrasts_postcue = log10(unique(data.contrast_postcue));
% no accuracy measure for orientation discrim when stim is absent = -1; no keyboard input = -2

idx.good = data.goodTrial==1;
for i = 1:numel(a.uniqueContrasts_postcue)
    % All trials 
    contrasts =         data.contrast_postcue; 
    idx.contrast = contrasts == a.uniqueContrasts_postcue(i);

    % contrastsPrecue = data.contrast_precue;
    % idx.contrastPrecue = contrastsPrecue == a.uniqueContrasts_postcue(i);
    % idxFPrecue = idx.contrastPrecue & idx.good; 
    
    if i==1 % contrast 0 
        idxF = idx.contrast & idx.good & idx.stimAbsent; 
    else
        idxF = idx.contrast & idx.good & idx.stimPresent; 
    end
    a.correctDis_allTrials{i}     = data.correctDis(idxF);
    a.correctDet_allTrials{i}     = data.correctDet(idxF);
    a.cueValidity_allTrials{i}    = data.cueValidity(idxF);
    
    a.sawGrating_allTrials{i}     = idx.sawGrating(idxF); 
    a.sawGratingNO_allTrials{i}   = idx.sawGratingNO(idxF); 
    a.sawOri_allTrials{i}         = idx.sawOri(idxF); 
    
    % Trials by attention 
    for iAtt = 1:numel(a.attHeaders)
        idx.att = data.cueValidity==a.attConds(iAtt);
        if i==1 % contrast 0
            idxF = idx.contrast & idx.good & idx.stimAbsent & idx.att;
        else
            idxF = idx.contrast & idx.good & idx.stimPresent & idx.att;
        end
        a.stimPresent{i,iAtt}     = idx.stimPresent(idxF); 
        a.goodTrial{i,iAtt}       = data.goodTrial(idxF);
        a.cueValidity{i,iAtt}     = data.cueValidity(idxF);

        a.correctDis{i,iAtt}      = data.correctDis(idxF);
        a.correctDet{i,iAtt}      = data.correctDet(idxF);
          
        % Detection 
        a.sawGrating{i,iAtt} = idx.sawGrating(idx.contrast & idx.good & idx.att);
        a.sawGratingNO{i,iAtt}    = idx.sawGratingNO(idxF); 
        a.sawOri{i,iAtt}          = idx.sawOri(idxF); 

        a.sawGrating_Present{i,iAtt}      = idx.sawGrating(idx.contrast & idx.good & idx.stimPresent & idx.att); 
        a.sawGrating_Absent{i,iAtt}       = idx.sawGrating(idx.good & idx.stimAbsent & idx.att); 

        a.sawOri_Present{i,iAtt}      = idx.sawOri(idx.contrast & idx.good & idx.stimPresent & idx.att); 
        a.sawOri_Absent{i,iAtt}       = idx.sawOri(idx.good & idx.stimAbsent & idx.att); 

        a.sawGratingNO_Present{i,iAtt}      = idx.sawGratingNO(idx.contrast & idx.good & idx.stimPresent & idx.att); 
        a.sawGratingNO_Absent{i,iAtt}       = idx.sawGratingNO(idx.good & idx.stimAbsent & idx.att); 
        
        % Discrimination 
        a.stimOri{i,iAtt}         = data.stimID_postcue(idxF); % 0-->absent, 1-->vertical, 2-->horizontal
        a.responseOri{i,iAtt}     = data.respDis(idxF); % 123-->H, 456-->V
        
        % Detection (all trials figure present and absent) to calculate d'
        % and c measures 
        a.sawGrating_absent{i,iAtt}       = idx.sawGrating_absent(idx.contrast & idx.good & idx.att & idx.stimAbsent); 
        a.sawGratingNO_absent{i,iAtt}     = idx.sawGratingNO_absent(idx.contrast & idx.good & idx.att & idx.stimAbsent); 
        a.sawOri_absent{i,iAtt}           = idx.sawOri_absent(idx.contrast & idx.good & idx.att & idx.stimAbsent); 

        % Detection (both figure present and absent trials) 
        % to calculate d' and c measures 
        a.sawGrating_presentAbsent{i,iAtt}       = idx.sawGrating(idx.contrast & idx.good & idx.att); 
        a.sawGratingNO_presentAbsent{i,iAtt}     = idx.sawGratingNO(idx.contrast & idx.good & idx.att); 
        a.sawOri_presentAbsent{i,iAtt}           = idx.sawOri(idx.contrast & idx.good & idx.att); 
        
        % For calcluating discrimination SDT measures 
        a.dis_correct_V{i,iAtt} = data.correctDis(idx.contrast & idx.att & idx.good & idx.stimPresent & data.stimID_postcue==1); 
        a.dis_correct_H{i,iAtt} = data.correctDis(idx.contrast & idx.att & idx.good & idx.stimPresent & data.stimID_postcue==2); 
    
        % Conditioning responses on the precue 
        a.correctDet_conditionedOnPrecue{i,iAtt} = data.correctDet_conditionedOnPrecue(idxF); 
        a.correctDis_conditionedOnPrecue{i,iAtt} = data.correctDis_conditionedOnPrecue(idxF); 

        % Retrieve the distractor contrasts 
        a.contrast_distractors{i,iAtt} = data.contrast_distractors(:,idxF);
        a.contrast_distractorsMax{i,iAtt} = max(data.contrast_distractors(:,idxF)); 
    
        % RT 
        a.RT{i,iAtt} = data.RT(idxF); 

        % Responses conditioned on seeing target
        a.correctDis_sawGrating{i,iAtt} = data.correctDis(idxF & idx.sawGrating);
        a.correctDis_sawOri{i,iAtt} = data.correctDis(idxF & idx.sawOri); 

        a.correctDet_sawGrating{i,iAtt} = data.correctDet(idxF & idx.sawGrating);
        a.correctDet_sawOri{i,iAtt} = data.correctDet(idxF & idx.sawOri); 
    end
end

for i = 1:numel(cFieldNames)
    n.(cVarNames{i}) = sum(cellfun(@sum, a.(cVarNames{i})),'all','omitnan'); % get detection resp type counts 
end

%% Calculate SDT measures
clear sdt
loglinearcorrection = 1; % to correct for floor or ceiling cells

% Detection SDT --> collapsed across for absent 
for iC = 1:numel(cFieldNames)
    for iL = 2:numel(a.uniqueContrasts_postcue)
        for iAtt = 1:numel(a.attHeaders)
            condFA = []; nfa = []; nnoise = [];
            condH = []; nh = []; nsignal = [];
            
            condName = sprintf('%s_Absent',cVarNames{iC}); 
            % condName = 'sawGrating_Absent';
            condFA = a.(condName){iL,iAtt}; % same across all stim strengths, this is storing reports of seeing when the contrast is 0
            nfa = numel(find(condFA==1));
            nnoise = numel(condFA); 
            
            condName = sprintf('%s_Present',cVarNames{iC}); 
            % condName = 'sawGrating_Present';
            condH = a.(condName){iL,iAtt};
            nh = numel(find(condH==1));
            nsignal = numel(condH);
            
            condName = sprintf('%s',cVarNames{iC}); 
            [sdt.(condName).dprimeDetect(iL,iAtt), sdt.(condName).criterionDetect(iL,iAtt)] = kt_dprime(nh,nfa,nsignal,nnoise,loglinearcorrection);
            
            sdt.(condName).det_nh(iL,iAtt) = nh;
            sdt.(condName).det_nfa(iL,iAtt) = nfa;
            sdt.(condName).det_nsignal(iL,iAtt) = nsignal;
            sdt.(condName).det_nnoise(iL,iAtt) = nnoise;
            sdt.(condName).det_ph(iL,iAtt) = nh/nsignal;
            sdt.(condName).det_pfa(iL,iAtt) = nfa/nnoise;
        end
    end
end

% Detection SDT (UV) - treats the saw grating and saw feature as graded
% responses along a common decision axis 

% Discrimination SDT 
for iL = 1:numel(a.uniqueContrasts_postcue)
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
var = a.sawGrating; 
varName = 'C1'; % 
varShortName = 'saw grating';
taskType = 'detection';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Detection, C2, p(saw shape)
p.iPlot = p.iPlot+1; 
var = a.sawOri; 
varName = 'C2'; 
varShortName = 'saw ori';
taskType = 'detection';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

% annotation = sprintf('%d of %d trials (%0.2f)',n.(cVarNames{2}),n.(cVarNames{1}),n.(cVarNames{2})/n.(cVarNames{1}));
% t = title(annotation);
% t.FontWeight = 'normal'; 

%% Detection, CX, p(saw figure and not saw shape)
p.iPlot = p.iPlot+1; 
var = a.sawGratingNO; 
varName = 'CX'; 
varShortName = 'saw grating and not ori';
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
if plotPSI
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
end

%% Plots AUC 
if plotPSI
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
end

%% Save analyses 
% General data parameters 
a.exptName         = s.exptName; 
a.exptShortName    = s.exptShortName;
a.site             = s.site; 
a.exptFolder       = s.exptFolder; 
a.gratingType      = s.gratingType; 
a.dataFolder       = s.dataFolder; 
a.subjectFolder    = s.subjectFolder; 
a.subjectIdx       = subjectIdx; 
a.subjectData      = data; 


