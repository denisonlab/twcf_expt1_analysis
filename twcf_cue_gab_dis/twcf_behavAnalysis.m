function [a,PSI] = twcf_behavAnalysis(dataAll, subjectIDs, subjectID, s, plotPSI)
% [a] = behav_analysis(dataAll,subjectIDs,subjectID)

% 1.4 cue gab dis  
% runBehav
% Inputs 
%   dataAll is structure of concatenated session data 
%   subjectIDs is cell array of subject IDs 
%   subjectID is string of subject of interest 
%   s is setup structure
% Outputs
%   a is analysis structure 

%% Do analysis on figure present or figure absent trials (move to setup) 
s.figureType     = 'figure present';  % figure means grating here 

%% Setup 
addpath(genpath(pwd))

% Get figure and psychometric function fitting parameters 
p = twcf_analysisParams; 

%% Load data 
subjectIdx = find(ismember(subjectIDs, subjectID));
clear data 
data = dataAll(subjectIdx);
s.subjectFolder = sprintf('%s/%s/%s/%s',s.dataFolder,s.exptShortName,subjectID,s.exptFolder);

% Make figure directory 
if s.saveFigs % change when using UCI+BU data to unlocked folder 
    figFolder = sprintf('%s/figs',s.subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% Set different conditions of interest and condition names 
% C
cVarNames      = {'testStronger'}; % 'sawGrating','sawOri','sawGratingNO'
cFieldNames    = {'C'}; 
cVarLabelNames = {'p(test stronger)'}; 

% A 
aVarNames      = {'correctDis'};
aFieldNames    = {'A'}; 
aVarLabelNames = {'p(correct discrimination)'}; 

% All var names 
varNames       = {'A','C'};

% Attention 
a.attConds   = unique(data.cueValidity); % -1 --> invalid, 1 --> valid, 0 --> neutral
a.attHeaders = {'invalid','neutral','valid'};

%% Get indices stim present, saw figure, saw shape, saw figure | saw shape from response buttons 
idx.stimPresent = logical(data.stimID_postcue(:)>0)'; % all trials present 1 or 2 only 

% data.resp is not coding the key press... is coding category of resp (see
% param table) 
% idx.sawGrating           = double(logical(data.resp==1 | data.resp==2 | data.resp==5 | data.resp==6)); % if saw shape, also saw figure 
% idx.sawGratingNO         = double(logical(data.resp==2 | data.resp==5)); % saw figure and not shape 
% idx.ratingDis          = double(logical(data.ratingDis); % responded test stronger 
idx.testStronger = data.ratingDis==1; 
idx.sawTilt = data.resp2==1; 

switch s.figureType 
    case 'figure present'
        idxCase = ~idx.stimPresent; 
    case 'figure absent'
        idxCase = idx.stimPresent; 
end
% idx.sawGrating(idxCase)           = NaN;
% idx.sawGratingNO(idxCase)         = NaN;
idx.sawOri(idxCase)               = NaN; % NaN orientation response when target absent

%% Break down data by line length and attention conditions 
a.contrasts_postcue = unique(data.contrast_postcue); % 7 stimulus strengths 
a.log10_contrasts_postcue = log10(unique(data.contrast_postcue));
% no accuracy measure for orientation discrim when stim is absent = -1; no keyboard input = -2

idx.good = data.goodTrial==1;
for i = 1:numel(a.contrasts_postcue)
    % all trials 
    contrasts     =         data.contrast_postcue; 
    idx.contrast = contrasts == a.contrasts_postcue(i);

    a.correctDis_allTrials{i}     = data.correctDis(idx.contrast & idx.good & data.correctDis~=-1 & idx.stimPresent);
    a.ratingDis_allTrials{i}      = data.ratingDis(idx.contrast & idx.good & idx.stimPresent);
    a.cueValidity_allTrials{i}    = data.cueValidity(idx.contrast & idx.good & idx.stimPresent);
    
    a.testStronger_allTrials{i}      = idx.testStronger(idx.contrast & idx.good & idx.stimPresent);  
    
    % trials by attention 
    for iAtt = 1:numel(a.attHeaders)
        idx.att = data.cueValidity==a.attConds(iAtt);
        a.stimPresent{i,iAtt}     = idx.stimPresent(idx.contrast & idx.good & idx.att & idx.stimPresent); 
        a.goodTrial{i,iAtt}       = data.goodTrial(idx.contrast & idx.att & idx.stimPresent);
        a.cueValidity{i,iAtt}     = data.cueValidity(idx.contrast & idx.good & idx.att & idx.stimPresent);

        a.correctDis{i,iAtt}      = data.correctDis(idx.contrast & idx.good & data.correctDis~=-1 & idx.att & idx.stimPresent);
        a.ratingDis{i,iAtt}       = data.ratingDis(idx.contrast & idx.good & idx.att & idx.stimPresent);
        
        % Subjective (strength) 
        a.testStronger{i,iAtt}        = idx.testStronger(idx.contrast & idx.good & idx.att & idx.stimPresent); 
        
        if i<4
            a.correctComparison{i,iAtt}   = ~idx.testStronger(idx.contrast & idx.good & idx.att & idx.stimPresent);
        elseif i>4
            a.correctComparison{i,iAtt}   = idx.testStronger(idx.contrast & idx.good & idx.att & idx.stimPresent); 
        end
           
        % Second repsonse (saw tilt or guess tilt) 
        a.sawTilt{i,iAtt}             = data.resp2(idx.contrast & idx.good & idx.att & idx.stimPresent); 
        % Test stronger | saw tilt 
        a.testStrongerSawTilt{i,iAtt} = idx.testStronger(idx.sawTilt & idx.contrast & idx.good & idx.att & idx.stimPresent); 
        
        % Objective (orientation discrimination) 
        a.stimOri{i,iAtt}        = data.stimID_postcue(idx.contrast & idx.att & idx.good & idx.stimPresent); % 0-->absent, 1-->vertical, 2-->horizontal
        a.responseOri{i,iAtt}    = data.respDis(idx.contrast & idx.att &idx.good & idx.stimPresent); % 123-->H, 456-->V  
        
        % For calcluating discrimination SDT measures 
        a.dis_correct_V{i,iAtt} = data.correctDis(idx.contrast & idx.att & idx.good & idx.stimPresent & data.stimID_postcue==1); 
        a.dis_correct_H{i,iAtt} = data.correctDis(idx.contrast & idx.att & idx.good & idx.stimPresent & data.stimID_postcue==2); 

        % RT 
        a.RT{i,iAtt}             = data.RT(idx.contrast & idx.att &idx.good & idx.stimPresent); 
    end
end

for i = 1:numel(cFieldNames)
    n.(cVarNames{i}) = sum(cellfun(@sum, a.(cVarNames{i})),'all','omitnan'); % get detection resp type counts 
end

%% Calculate SDT measures
clear sdt

% Detection SDT: taking true reference strength to categorize H FA
% Ignore when stimuli are equal? 
% recalculate collapsed across line lengths for absent 

condName = {'comparison','comparisonByLL'}; % comparison by line length
for iC = 1:numel(condName)
    for iAtt = 1:numel(a.attHeaders)
        condFA = []; nfa = []; nnoise = [];
        condH = []; nh = []; nsignal = [];
    
        switch condName{iC}
            case 'comparison'
                for iL = 1:3 % False alarms; reporting test stronger when it was weaker than reference
                    condFA = cat(2,condFA,a.testStronger{iL,iAtt});
                end
                nfa = sum(condFA);
                nnoise = numel(condFA);

                for iL = 5:7 % Hits; reporting test stronger when it was stronger than reference
                    condH = cat(2,condH,a.testStronger{iL,iAtt});
                end
                nh = sum(condH);
                nsignal = numel(condH);
                [sdt.(condName{iC}).dprimeDetect(iAtt), sdt.(condName{iC}).criterionDetect(iAtt)] = rd_dprime2(nh,nfa,nsignal,nnoise);

                sdt.(condName{iC}).det_nh(iAtt) = nh;
                sdt.(condName{iC}).det_nfa(iAtt) = nfa;
                sdt.(condName{iC}).det_nsignal(iAtt) = nsignal;
                sdt.(condName{iC}).det_nnoise(iAtt) = nnoise;
                sdt.(condName{iC}).det_ph(iAtt) = nh/nsignal;
                sdt.(condName{iC}).det_pfa(iAtt) = nfa/nnoise;

            case 'comparisonByLL'
                for iL = 1:7
                    if iL==4
                        sdt.(condName{iC}).dprimeDetect(iL,iAtt) = NaN;
                        sdt.(condName{iC}).criterionDetect(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_nh(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_nfa(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_nsignal(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_nnoise(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_ph(iL,iAtt) = NaN;
                        sdt.(condName{iC}).det_pfa(iL,iAtt) = NaN;
                    else
                        if iL<4
                            condFA = a.testStronger{iL,iAtt};

                            nfa = sum(condFA);
                            nnoise = numel(condFA); 

                            % loglinear adjustment (Stanislaw & Todorov 1999) 
                            nh = 0.5;
                            nsignal = 1; 

                        elseif iL>4
                            condH = a.testStronger{iL,iAtt};

                            nfa = 0.5; 
                            nnoise = 1; 

                            nh = sum(condH);
                            nsignal = numel(condH); 
                        end
                        % nfa = sum(condFA);
                        % nnoise = numel(condFA);
                        % nh = sum(condH);
                        % nsignal = numel(condH);
                        loglinearcorrection = 1; 
                        [sdt.(condName{iC}).dprimeDetect(iL,iAtt), sdt.(condName{iC}).criterionDetect(iL,iAtt)] = kt_dprime(nh,nfa,nsignal,nnoise,loglinearcorrection);

                        sdt.(condName{iC}).det_nh(iL,iAtt) = nh;
                        sdt.(condName{iC}).det_nfa(iL,iAtt) = nfa;
                        sdt.(condName{iC}).det_nsignal(iL,iAtt) = nsignal;
                        sdt.(condName{iC}).det_nnoise(iL,iAtt) = nnoise;
                        sdt.(condName{iC}).det_ph(iL,iAtt) = nh/nsignal;
                        sdt.(condName{iC}).det_pfa(iL,iAtt) = nfa/nnoise;
                    end
                end
        end
    end
end

% Discrimination sdt 
for iL = 1:numel(a.contrasts_postcue)
    for iAtt = 1:numel(a.attHeaders)
        condFA = []; nfa = []; nnoise = []; 
        condH = []; nh = []; nsignal = []; 

        condName = 'dis_correct_V'; % take one orientation to be signal, other as noise 

        condFA = a.(condName){iL,iAtt}; 
        nfa = numel(find(condFA==0)); 
        nnoise = numel(condFA); 

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
sgtitle(sprintf('expt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.sites{subjectIdx}, subjectID, max(data.sessions)))
p.iPlot = 0;

%% Discrimination, A, p(correct)
p.iPlot = p.iPlot+1; 
var = a.correctDis; 
varName = 'A'; % discriminationpCorrect
varShortName = 'correct';
taskType = 'discrimination';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Subjective rating, C, p(test stronger)
p.iPlot = p.iPlot+1; 
var = a.testStronger; 
varName = 'C1'; 
varShortName = 'test stronger';
taskType = 'subjective';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Subjective rating, C, p(saw tilt)
p.iPlot = p.iPlot+1; 
var = a.sawTilt; 
varName = 'C2'; 
varShortName = 'saw tilt';
taskType = 'subjective';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Subjective rating, C, p(test stronger | saw tilt)
p.iPlot = p.iPlot+1; 
var = a.testStrongerSawTilt; 
varName = 'C3'; 
varShortName = 'test stronger | saw tilt';
taskType = 'subjective';

a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Detection, C1, p(saw figure)
% p.iPlot = p.iPlot + 1; 
% var = a.sawGrating; 
% varName = 'C1'; % 
% varShortName = 'saw grating';
% taskType = 'detection';
% 
% a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

%% Detection, C2, p(saw shape)
% p.iPlot = p.iPlot+1; 
% var = a.sawOri; 
% varName = 'C2'; 
% varShortName = 'saw ori';
% taskType = 'detection';
% 
% a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

% annotation = sprintf('%d of %d trials (%0.2f)',n.(cVarNames{2}),n.(cVarNames{1}),n.(cVarNames{2})/n.(cVarNames{1}));
% t = title(annotation);
% t.FontWeight = 'normal'; 

%% Detection, CX, p(saw figure and not saw shape)
% p.iPlot = p.iPlot+1; 
% var = a.sawGratingNO; 
% varName = 'CX'; 
% varShortName = 'saw grating and not ori';
% taskType = 'detection';
% 
% a = twcf_plotFit(var,varName,varShortName,a,p,taskType,s);

% annotation = sprintf('%d of %d trials (%0.2f)',n.(cVarNames{3}),n.(cVarNames{1}),n.(cVarNames{3})/n.(cVarNames{1}));
% t = title(annotation);
% t.FontWeight = 'normal'; 

%% Save figure with detection discrimination subplots 
if s.saveFigs 
    figTitle = sprintf('%s_%s_%s_%s_CAbyXAtt',s.exptFolder,s.sites{subjectIdx},subjectID,datestr(now,'yymmdd')); 
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType)) 
end

%% Plots fitted params (subject-level)
varNames       = {'A','C1'};
figure
set(gcf,'Position',[0 0 400 p.style.figHeight/4 * numel(varNames)])
sgtitle(sprintf('Fitted params\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.sites{subjectIdx}, subjectID, max(data.sessions)))
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
    figTitle = sprintf('%s_%s_%s_%s_fittedParams',s.exptFolder,s.sites{subjectIdx},subjectID,datestr(now,'yymmdd'));
    saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType))
end

%% Plots CvA
% Testing performance-relative subjective inflation (PSI)
cFieldNames    = {'C1','C2','C3'}; 
if plotPSI
    figure
    set(gcf,'Position',[0 0 350 p.style.figHeight])
    sgtitle(sprintf('Testing PSI - Metaperceptual function (CvA)\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.sites{subjectIdx}, subjectID, max(data.sessions)))
    p.iPlot = 1;

    varA = a.A;
    for iC = 1:numel(cFieldNames)
        varC = a.(cFieldNames{iC});
        PSI.(cFieldNames{iC}) = twcf_plotCvA(varA,varC,p,a);
        p.iPlot = p.iPlot+1;
    end

    if s.saveFigs
        figTitle = sprintf('%s_%s_%s_%s_CvA',s.exptFolder,s.sites{subjectIdx},subjectID,datestr(now,'yymmdd'));
        saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType))
    end
end

%% Plots AUC and relative AUC
% constrain by maximum possible area 
% maximum shared x range * 1 
if plotPSI
    figure
    set(gcf,'Position',[0 0 350 (p.style.figHeight/3)*2])
    sgtitle(sprintf('Testing PSI - relative AUC of CvA\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.sites{subjectIdx}, subjectID, max(data.sessions)))
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
        figTitle = sprintf('%s_%s_%s_%s_AUC',s.exptFolder,s.sites{subjectIdx},subjectID,datestr(now,'yymmdd'));
        saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType ))
    end
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



