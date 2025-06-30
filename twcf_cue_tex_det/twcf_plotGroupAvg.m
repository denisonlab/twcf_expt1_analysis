function [g2] = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g)

% 1.1 cue tex detection 
% var is structure of group analyses (from behav_analysis) 

%% Settings 
annotateStats = 1; 
annotateN = 1; 

%% Setup 
p = twcf_analysisParams;
set(0,'DefaultAxesTitleFontWeight','normal');

% Figure directory 
% figDir = sprintf('%s/%s', s.groupFigDir, datestr(now,'yymmdd')); 
figDir = sprintf('%s/%s', s.groupFigDir, 'manuscript_pdf'); 
% figDir = sprintf('%s/%s', s.groupFigDir, 'analysis'); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

plotFigAbsent = 0; % plots FAR

%% Plot labels 
cVarNames = {'C1','C2'}; % ,'CX','CY'}; 
% cVarShortNames = {'saw figure','saw shape','saw figure not shape','saw shape only'}; 
% cVarLabelNames = {'p("saw figure")','p("saw shape")','p("saw figure not shape")','p("saw shape only")'};
cVarShortNames = {'saw stimulus','saw feature'}; 
cVarLabelNames = {'p("saw stimulus")','p("saw feature")'};

% A 
aVarNames      = {'correctDis'};
aFieldNames    = {'A'}; 
aVarLabelNames = {'p(correct discrimination)'}; 

% clear a.attHeaders; 
% a.attConds   = {-1,0,1};  % -1 --> invalid, 1 --> valid, 0 --> neutral
% a.attHeaders = {'invalid','neutral','valid'};

%% Average group data
val_x = []; val_logx = []; val_C = []; val_A = []; 
val_stimPresent = []; 

for iS = 1:numel(subjectIDs)
    val_x(iS,:)      = a(iS).uniqueLineLengths_inDeg_postcue;
    val_logx(iS,:)   = a(iS).log10_uniqueLineLengths_inDeg_postcue;
    val_C(iS,:,:)    = a(iS).(cFieldNames{iC}).yMean;
    val_A(iS,:,:)    = a(iS).A.yMean;
    for iL = 1:7
        for iAtt = 1:3
            val_stimPresent(iS,iL,iAtt) = sum(a(iS).stimPresent{iL,iAtt});
            val_RT(iS,iL,iAtt)   = mean(a(iS).RT{iL,iAtt},'omitnan'); 
        end
    end
end

val_C_absent = g(iC).det_pfa.val;
val_C_dprime = g(iC).dprimeDetect.val; 
val_C_criterion = g(iC).criterionDetect.val; 

%% Collect group data 
group_a.attConds = a(1).attConds; 
group_a.attHeaders = a(1).attHeaders; 
group_a.uniqueLineLengths_inDeg_postcue = mean(val_x,1); 
group_a.uniqueLineLengths_inDeg_postcue_S = val_x; 
group_a.log10_uniqueLineLengths_inDeg_postcue = mean(val_logx,1); 
group_a.outOfNum = squeeze(max(val_stimPresent,[],1)); 

%% Discrimination, A, p(correct)
figure
set(gcf,'Position',[0 0 380 280])
% sgtitle(sprintf('expt: %s | site: %s', s.exptName, s.site))
p.iPlot = 0; 

p.iPlot = p.iPlot+1;
s.figureType = 'figure present'; 
varSubjects = val_A; % subjects x strength x att 
var = squeeze(mean(val_A,1));  
varName = 'A'; % discriminationpCorrect
varShortName = 'correct discrimination';
taskType = 'discrimination';

a_A = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,varSubjects);
ax = gca; 
sizeAx = plotboxpos(ax);

% --- Add n annotation ---
if annotateN
    twcf_annotateN(a)
end

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_A(x)',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end

%% Detection, C1 or 2, p(saw shape), p(saw figure) 
figure
set(gcf,'Position',[p.size.rect])
% sgtitle(sprintf('expt: %s | site: %s', s.exptName, s.site))
p.iPlot = 0; 

p.iPlot = p.iPlot+1; 
var = squeeze(mean(val_C,1)); 
varName = cVarNames{iC}; 
varShortName = cVarShortNames{iC};
taskType = 'detection';

a_C.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C);

if plotFigAbsent
    var = squeeze(mean(val_C_absent,1)); % average subjects 
    varName = cVarNames{iC};
    varShortName = cVarShortNames{iC};
    taskType = 'detection';
    a_C_absent.(cVarNames{iC}) = twcf_plotFitGroupAbsent(var,varName,varShortName,group_a,p,taskType,s);
end

if s.saveFigs
    if plotFigAbsent
        figTitle = sprintf('group_n%d_%s_C%d(x)_absent',numel(subjectIDs),s.site,iC);
    else
        figTitle = sprintf('group_n%d_%s_C%d(x)',numel(subjectIDs),s.site,iC);
    end
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end

%% Detection d' (overlay stimulus and feature detection) 
if iC==1
    % Stimulus detection
    val_C_dprime_stimulus = g(1).dprimeDetect.val; 
    val_C_criterion_stimulus = g(1).criterionDetect.val;
    % Feature detection 
    val_C_dprime_feature = g(2).dprimeDetect.val; 
    val_C_criterion_feature = g(2).criterionDetect.val;

    figure
    set(gcf,'Position',p.size.rect)
    p.iPlot = 0;

    p.iPlot = p.iPlot+1;

    % plot stimulus 
    var = squeeze(mean(val_C_dprime_stimulus,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
    taskType = 'dprime_stimulus';
    twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);

    % plot feature 
    var = squeeze(mean(val_C_dprime_feature,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
    taskType = 'dprime_feature';
    twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
    ylabel({'Detection \itd''\rm'},'FontSize',p.style.textAxisSize)
    ylim([-0.5 4])

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_dprime_overlaid',numel(subjectIDs),s.site);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Detection c (overlay stimulus and feature detection) 
if iC==1
    figure
    set(gcf,'Position',p.size.rect)
    p.iPlot = 0;

    p.iPlot = p.iPlot+1;

    % plot stimulus 
    var = squeeze(mean(val_C_criterion_stimulus,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
    taskType = 'criterion_stimulus';
    twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);

    % plot feature 
    var = squeeze(mean(val_C_criterion_feature,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
    taskType = 'criterion_feature';
    twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
    ylabel({'Detection criterion'},'FontSize',p.style.textAxisSize)
    yticks([-1:1:2])
    yline(0)

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_criterion_overlaid',numel(subjectIDs),s.site);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end


%% Detection d' by detection criterion 
if iC==1||2
    figure
    figureStyle 
    set(gcf,'Position',[0 0 380 280])

    varX = squeeze(mean(val_C_dprime,1));
    varY = squeeze(mean(val_C_criterion,1));

    xContrasts = 1:7;

    % Lines first 
    for iAtt = 1:numel(a(1).attHeaders)
        plot(varX(xContrasts,iAtt),varY(xContrasts,iAtt),'Color',p.style.attColors(iAtt,:),'linewidth',p.style.fitLineWidth)
    end

    % Scatter second 
    for iAtt = 1:numel(a(1).attHeaders)
        x = varX(xContrasts,iAtt);
        y = varY(xContrasts,iAtt);

        color = p.style.attColors(iAtt,:); 
        sc = scatter(x,y,p.style.sz,'MarkerFaceColor',color,'MarkerFaceAlpha',1,'MarkerEdgeColor','w','LineWidth',2);
        sc.Marker = 'o';
    end

    % a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_criterion);
    ylabel({'Detection criterion'},'FontSize',p.style.textAxisSize)
    xlabel('Detection {\itd''}')

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_criterion_dprime',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Discrimination d' (conditioned on when the stimulus was present)
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0;

p.iPlot = p.iPlot+1;
var = squeeze(mean(g(1).dprimeDiscrim.val,1));
varName = 'discrimdprime';
varShortName = sprintf('Discrimination d''\n%s');
taskType = 'dprime';

a_A_dprime = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,g(1).dprimeDiscrim.val);
ylabel({'Discrimination \itd''\rm'},'FontSize',p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_discrimination_dprime',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end

%% Discrimination criterion
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0;

p.iPlot = p.iPlot+1;
var = squeeze(mean(abs(g(1).criterionDiscrim.val),1));
varName = 'discrimcriterion';
varShortName = sprintf('Discrimination criterion\n%s',cVarShortNames{iC});
taskType = 'criterion';

a_A_criterion.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,abs(val_C_criterion));
ylabel({'Abs (Discrimination criterion)'},'FontSize',p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_discrimination_criterion',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end


%% D' UV/EV SDT
% load gsdt structure loaded
load('1.1sdt.mat')
varName = cVarNames{iC};
group_a.outOfNum = squeeze(a(1).(varName).OutOfNum)';
s.figureType = 'figure present';

sdtType = 'UV'; % EV, UV
switch sdtType 
    case 'EV' % equal variance
        % Stimulus detection
        val_C_dprime_stimulus = gsdt.dprime(:,:,:,1);
        val_C_criterion_stimulus = gsdt.criterion(:,:,:,1);
        % Feature detection
        val_C_dprime_feature = gsdt.dprime(:,:,:,2);
        val_C_criterion_feature = gsdt.criterion(:,:,:,2);
    case 'UV'
        % Stimulus detection
        val_C_dprime_stimulus = gsdt.d_a(:,:,:,1);
        val_C_criterion_stimulus = gsdt.c_a(:,:,:,1);
        % Feature detection
        val_C_dprime_feature = gsdt.d_a(:,:,:,2);
        val_C_criterion_feature = gsdt.c_a(:,:,:,2);
end

figure
set(gcf,'Position',p.size.rect)
p.iPlot = 0;

p.iPlot = p.iPlot+1;

% plot stimulus
var = squeeze(mean(val_C_dprime_stimulus,1,'omitnan'));
varName = cVarNames{iC};
cVarShortNames = {'saw figure','saw orientation'};
varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
switch sdtType 
    case 'EV'
        taskType = 'dprime_stimulus';
    case 'UV'
        taskType = 'dprime_stimulus_UV';
end
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);

% plot feature
var = squeeze(mean(val_C_dprime_feature,1,'omitnan'));
varName = cVarNames{iC};
cVarShortNames = {'saw grating','saw orientation'};
varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
switch sdtType
    case 'EV'
        taskType = 'dprime_feature';
    case 'UV'
        taskType = 'dprime_feature_UV';
end
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
switch sdtType
    case 'EV'
        ylabel({'Detection \itd''\rm'},'FontSize',p.style.textAxisSize)
    case 'UV'
        ylabel({'\itd_{a}'},'FontSize',p.style.textAxisSize)
end
ylim([-0.5 4])
xlim([0.0694 0.68])

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_dprime_overlaid_%s',numel(subjectIDs),s.site,sdtType);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end

%% Criterion UV/EV SDT 
figure
set(gcf,'Position',p.size.rect)
p.iPlot = 0;

p.iPlot = p.iPlot+1;

% plot stimulus
var = squeeze(mean(val_C_criterion,1,'omitnan'));
varName = cVarNames{iC};
cVarShortNames = {'saw grating','saw orientation'};
varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
switch sdtType 
    case 'EV'
        taskType = 'criterion_stimulus';
    case 'UV'
        taskType = 'criterion_stimulus_UV';
end
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);

% plot feature
var = squeeze(mean(val_C_criterion_feature,1,'omitnan'));
varName = cVarNames{iC};
cVarShortNames = {'saw grating','saw orientation'};
varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
switch sdtType 
    case 'EV'
        taskType = 'criterion_feature';    
    case 'UV'
        taskType = 'criterion_feature_UV';       
end
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
yticks([-1:1:2])
ylim([-1 2.5])
yline(0,'Color',[0.7 0.7 0.7])
xlim([0.0694 0.68])

switch sdtType
    case 'EV'
        ylabel('Detection criterion','FontSize',p.style.textAxisSize)
    case 'UV'
        ylabel('\itc_{a}','FontSize',p.style.textAxisSize)
end

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_criterion_overlaid_%s',numel(subjectIDs),s.site,sdtType);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end


%% Slope UV SDT 
switch sdtType
    case 'UV'
        figure
        set(gcf,'Position',p.size.rect)
        p.iPlot = 0;

        p.iPlot = p.iPlot+1;

        % plot slope
        var = squeeze(mean(gsdt.s,1,'omitnan'));
        taskType = 'slope_UV'; % borrow this setting
        twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);

        ylim([0.5 2])
        yticks([-2:0.5:2])
        yl = yline(1,'Color',[187 187 187]/255); 
        uistack(yl,'bottom'); % why doesn't this work
        ylabel('Slope','FontSize',p.style.textAxisSize)
        xlim([0.0694 0.68])
        
        if s.saveFigs
            figTitle = sprintf('group_n%d_%s_slope',numel(subjectIDs),s.site);
            export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
        end
end

%% Detection d' by detection criterion 
if iC==1||2
    figure
    figureStyle
    set(gcf,'Position',[0 0 370 280]) % p.size.rect

    for iV=1:2 % stimulus then feature
        clear varX varY se
        switch iV 
            case 1
                varX = squeeze(mean(val_C_dprime_stimulus,1));   % val_C_criterion_stimulus
                varY = squeeze(mean(val_C_criterion_stimulus,1));
                se = squeeze(std(val_C_criterion_stimulus,1))/sqrt(30);
            case 2
                varX = squeeze(mean(val_C_dprime_feature,1));   % val_C_criterion_stimulus
                varY = squeeze(mean(val_C_criterion_feature,1));
                se = squeeze(std(val_C_criterion_feature,1))/sqrt(30);
        end

        xContrasts = 1:7;
        % Lines first
        for iAtt = 1:numel(a(1).attHeaders)
            switch iV
                case 1
                    linestyle = ':';
                    color = p.style.attColorsMutedLight(iAtt,:);
                case 2
                    linestyle = '-';
                    color = p.style.attColors(iAtt,:);
            end
            
            errorbar(varX(xContrasts,iAtt),varY(xContrasts,iAtt),se(xContrasts,iAtt),'Color',color,'CapSize',0,'LineWidth',1,'LineStyle','none')
            plot(varX(xContrasts,iAtt),varY(xContrasts,iAtt),'Color',color,'linewidth',p.style.fitLineWidth,'linestyle',linestyle)
        end
    end

    ylabel({'Detection criterion'},'FontSize',p.style.textAxisSize)
    xlabel('Detection sensitivity','FontSize',p.style.textAxisSize)
    yticks([-1:1:2])
    yline(0,'Color',[0.7 0.7 0.7])

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_criterion_dprime',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Reaction time
figure
set(gcf,'Position',[0 0 380 280])

var = squeeze(mean(val_RT,1));  
varName = cVarNames{iC}; 
varShortName = 'RT';
taskType = 'RT';

group_a.outOfNum = a(1).(varName).OutOfNum'; % squeeze(max(val_stimPresent,[],1)); 

a_RT = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_RT,0,a);

% --- Add n annotation ---
if annotateN
    twcf_annotateN(a)
end

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_RT',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% Detection d'
if iC==1||2
    figure
    set(gcf,'Position',[0 0 380 280])
    p.iPlot = 0;

    p.iPlot = p.iPlot+1;
    var = squeeze(mean(val_C_dprime,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw figure','saw shape'};
    varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
    taskType = 'dprime';

    a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
    ylabel({'Detection \itd''\rm'},'FontSize',p.style.textAxisSize)

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_dprime',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Detection criterion
if iC==1||2
    figure
    set(gcf,'Position',[0 0 380 280])
    p.iPlot = 0;

    p.iPlot = p.iPlot+1;
    var = squeeze(mean(val_C_criterion,1));
    varName = cVarNames{iC};
    cVarShortNames = {'saw figure','saw shape'};
    varShortName = sprintf('Detection criterion\n%s',cVarShortNames{iC});
    taskType = 'criterion';

    a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_criterion);
    ylabel('Detection criterion','FontSize',p.style.textAxisSize)

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_criterion',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Detection criterion difference by attention 
if iC==1||2
    figure
    set(gcf,'Position',[0 0 380 280])
    p.iPlot = 0;

    p.iPlot = p.iPlot+1;
    var = squeeze(mean(val_C_criterion,1));
    varDiff_S = val_C_criterion(:,:,3) - val_C_criterion(:,:,1); % attended - unattended 
    varDiff = var(:,3)-var(:,1); % attended - unattended 
    varName = cVarNames{iC};
    % cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Criterion difference');
    taskType = 'criterion_diff';

    a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(varDiff,varName,varShortName,group_a,p,taskType,s,varDiff_S);
    % ylabel({'Detection criterion difference'},'FontSize',p.style.textAxisSize)

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_criterionDiff',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% Detection d' by detection criterion 
if iC==1||2
    figure
    figureStyle 
    set(gcf,'Position',[0 0 380 280])

    varX = squeeze(mean(val_C_dprime,1));
    varY = squeeze(mean(val_C_criterion,1));

    xContrasts = 1:7; % 2:8; 
    
    for iAtt = 1:numel(a(1).attHeaders)
        plot(varX(xContrasts,iAtt),varY(xContrasts,iAtt),'x')
    end

    xlabel('Detection {/itd''}')
    ylabel('Criterion')

    varName = cVarNames{iC};
    cVarShortNames = {'saw grating','saw orientation'};
    varShortName = sprintf('Detection criterion\n%s',cVarShortNames{iC});
    taskType = 'criterion';

    a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_criterion);
    ylabel({'Detection criterion'},'FontSize',p.style.textAxisSize)

    if s.saveFigs
        figTitle = sprintf('group_n%d_%s_C%d(x)_criterion_dprime',numel(subjectIDs),s.site,iC);
        export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
    end
end

%% CvA: p(correct) v p(seen)
figure
set(gcf,'Position',[0 0 380 280])
% sgtitle(sprintf('Testing PSI - Metaperceptual function (CvA)\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 1; 

varA = a_A.A;
varC = a_C.(cVarNames{iC}).(cVarNames{iC}); % squeeze(mean(val_C,1)); 
PSIgroup = twcf_plotCvAGroup(varA,varC,p,a);

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_C%dvA',numel(subjectIDs),s.site,iC);
    % saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10') 
end

%% plot AUC of CvA
figure
set(gcf,'Position',[0 0 200 280])
% sgtitle(sprintf('Testing PSI - Metaperceptual function (CvA)\nexpt: %s | site: %s \n subject: %s | %d sessions', s.exptName, s.site, subjectID, max(data.sessions)))
p.iPlot = 1; 

varA = a_A.A;
varC = a_C.(cVarNames{iC}).(cVarNames{iC}); 

subplot 111
subplot('Position', p.size.AUC)
hold on
figureStyle
for iAtt = 1:numel(varA.attHeaders)
    clear x; clear y;
    x = [1 2 3];
    % y = PSIgroup.AUC;
    y(iAtt) = mean(g(iC).AUC.val(:,iAtt),1,'omitnan'); 
    % scatter(x(iAtt),y(iAtt),p.style.sz,'filled','MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',1,'MarkerEdgeColor','w')
    bar(x(iAtt),y(iAtt),'FaceColor',p.style.attColorsMuted(iAtt,:),'EdgeColor','none','FaceAlpha',1)
    errb = errorbar( x(iAtt), y(iAtt), std(g(iC).AUC.val(:,iAtt))/sqrt(numel(g(iC).AUC.val(:,iAtt)))  );
    errb.LineWidth = 2; 
    errb.Color = 'k'; 
    errb.CapSize = 0; 
end
xlim([0 4])
% ylabel(sprintf('AUC of\n%s vs p(correct)',cVarLabelNames{iC}))
ylabel('AUC')
xticks([1,2,3])
xticklabels({'I','N','V'})

xAX = get(gca,'XAxis');  
yAX = get(gca,'YAxis');  
xlab = get(gca,'XLabel');
ylab = get(gca,'YLabel');

set(xAX,'FontSize', p.style.textAxisSize) % larger labels 
set(yAX,'FontSize', p.style.textTickSize)
set(xlab,'FontSize', p.style.textAxisSize)
set(ylab,'FontSize', p.style.textAxisSize)

% Stats annotation
if annotateStats
    if iC==1
        yPs = [0.88, 0.78, 0.68];
        ylim([0 0.32])
        yticks(0.:0.1:32)
        starText = {'***','***','**','*'}; % main effect, IV, IN, VN
    elseif iC==2
        yPs = [0.88, 0.78, 0.68];
        ylim([0 0.17])
        yticks(0.:0.05:17)
        starText = {'*','**','*'}; 
    elseif iC==3
        yPs = [0.88, 0.78, 0.68];
    end
    yl = ylim;
    xl = xlim;
    txt = twcf_annotateStats(2,max(yl)*p.sigStarMainEffect,starText{1}); % main effect
    
    if numel(starText)>1
        % I vs V
        yP = yPs(1);
        txt = twcf_annotateStats(2,max(yl)*yP,starText{2});
        twcf_drawBracket(1,3,yP)
    end

    if numel(starText)>2
        % I vs N
        yP = yPs(2);
        txt = twcf_annotateStats(1.5,max(yl)*yP,starText{3});
        twcf_drawBracket(1,2,yP)
    end

    if numel(starText)>3
        % N vs V
        yP = yPs(3);
        txt = twcf_annotateStats(2.5,max(yl)*yP,starText{4});
        twcf_drawBracket(2,3,yP)
    end
end

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_AUC_C%dvA',numel(subjectIDs),s.site,iC);
    % saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-dALLOWPSTRANSPARENCY','-dNOSAFER','-transparent','-p10') 
end

%% AUC of CvA by site 
figure
set(gcf,'Position',[0 0 330 280])

varA = a_A.A;
varC = a_C.(cVarNames{iC}).(cVarNames{iC}); 

subplot 111
subplot('Position', p.size.AUCbySite)
hold on
figureStyle
xPos = [1 2 3];
sitelabs = {'BU','UCI'}; 
for iSite = 1:2
    for iAtt = 1:numel(varA.attHeaders)
        clear x; clear y;
        if iSite==1 % BU 
            x = xPos; 
            yIdx = find(contains(sites,'BU')); 
        elseif iSite==2 % UCI 
            x = xPos+4;
            yIdx = find(contains(sites,'UCI')); 
        end
        y(iAtt) = mean(g(iC).AUC.val(yIdx,iAtt),1,'omitnan');
        bar(x(iAtt),y(iAtt),'FaceColor',p.style.attColorsMuted(iAtt,:),'EdgeColor','none','FaceAlpha',1)
        errb = errorbar( x(iAtt), y(iAtt), std(g(iC).AUC.val(yIdx,iAtt))/sqrt(numel(yIdx))  );
        errb.LineWidth = 2;
        errb.Color = 'k';
        errb.CapSize = 0;
    end
end
xlim([0 8])
% ylabel(sprintf('AUC of\n%s vs p(correct)',cVarLabelNames{iC}))
% yticks(0.:0.01:1)
ylabel('AUC')
ylim([0 0.32])

x = [1 2 3 5 6 7]; 
xticks(x)
xticklabels({'I','N','V','I','N','V'})

xAX = get(gca,'XAxis');  
yAX = get(gca,'YAxis');  
xlab = get(gca,'XLabel');
ylab = get(gca,'YLabel');

set(xAX,'FontSize', p.style.textAxisSize) % larger labels 
set(yAX,'FontSize', p.style.textTickSize)
set(xlab,'FontSize', p.style.textAxisSize)
set(ylab,'FontSize', p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_AUC_C%dvA_bySite',numel(subjectIDs),s.site,iC);
    % saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-dALLOWPSTRANSPARENCY','-dNOSAFER','-transparent','-p10') 
end

%% Save variables 
% g is n x val 
% g2 is data averaged to group x val 
g2.val_x    = val_x; 
g2.val_logx = val_logx; 
g2.val_C    = val_C; 
g2.val_A    = val_A; 

