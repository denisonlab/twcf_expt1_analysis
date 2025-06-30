function [g2] = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g)

% 1.4 cue gab discrimination 
% var is structure of group analyses (from behav_analysis) 

%% Settings 
annotateStats = 1; 
annotateN = 1; 

%% Setup 
p = twcf_analysisParams;
set(0,'DefaultAxesTitleFontWeight','normal');

% Figure directory 
% figDir = sprintf('%s/%s', s.groupFigDir, datestr(now,'yymmdd')); 
figDir = sprintf('%s/manuscript_pdf', s.groupFigDir); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

plotFigAbsent = 0; % Plots figure absent  

%% Plot labels 
cVarNames      = {'C1','C2','C3'}; 
cVarShortNames = {'test stronger','saw feature','test stronger & saw feature'}; 
cVarLabelNames = {'p("test stronger")','p("saw feature")','p("test stronger | saw feature")'};

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

nStrengths = 7; 
for iS = 1:numel(subjectIDs)
    val_x(iS,:)      = a(iS).contrasts_postcue;
    val_logx(iS,:)   = a(iS).log10_contrasts_postcue;
    val_C(iS,:,:)    = a(iS).(cFieldNames{iC}).yMean; % 7 x 3 
    val_A(iS,:,:)    = a(iS).A.yMean;
    for iL = 1:nStrengths
        for iAtt = 1:3
            val_stimPresent(iS,iL,iAtt) = sum(a(iS).stimPresent{iL,iAtt});
            val_RT(iS,iL,iAtt)   = mean(a(iS).RT{iL,iAtt},'omitnan'); 
        end
    end
end

for iS = 1:numel(subjectIDs)
    for iL = 1:nStrengths
        for iAtt = 1:3
            val_C_comparison(iS,iL,iAtt) = mean(PSI(iS).(cFieldNames{iC}).correctComparison{iL,iAtt},'omitnan');
        end
    end
end

% val_C_absent = g(iC).det_pfa.val;
% val_C_dprime = g(iC).dprimeDetect.val; 
% val_C_criterion = g(iC).criterionDetect.val; 

%% Collect group data 
group_a.attConds = a(1).attConds; 
group_a.attHeaders = a(1).attHeaders; 
group_a.uniqueContrasts_postcue = mean(val_x,1); 
group_a.log10_uniqueContrasts_postcue = mean(val_logx,1); 

%% Discrimination, A, p(correct)
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0; 

p.iPlot = p.iPlot+1;
s.figureType = 'figure present'; 
var = squeeze(mean(val_A,1));  
varName = 'A'; % discriminationpCorrect
varShortName = 'correct discrimination';
taskType = 'discrimination';

group_a.outOfNum = a(1).(varName).OutOfNum'; % squeeze(max(val_stimPresent,[],1)); 

a_A = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_A);

% --- Add n annotation ---
if annotateN
    twcf_annotateN(a)
end

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_A(x)',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% Detection, C1 or 2, p(test stronger) 
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0; 

p.iPlot = p.iPlot+1; 
var = squeeze(mean(val_C,1)); 
varName = cVarNames{iC}; 
varShortName = cVarShortNames{iC};
taskType = 'detection';

group_a.outOfNum = a(1).(varName).OutOfNum';

referenceAnnotation = 1; 
a_C.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C,referenceAnnotation);

% if plotFigAbsent
%     var = squeeze(mean(val_C_absent,1)); % average subjects 
%     varName = cVarNames{iC};
%     varShortName = cVarShortNames{iC};
%     taskType = 'detection';
%     a_C_absent.(cVarNames{iC}) = twcf_plotFitGroupAbsent(var,varName,varShortName,group_a,p,taskType,s);
% end

if s.saveFigs
    if plotFigAbsent
        figTitle = sprisntf('group_n%d_%s_C%d(x)_absent',numel(subjectIDs),s.site,iC);
    else
        figTitle = sprintf('group_n%d_%s_C%d(x)',numel(subjectIDs),s.site,iC);
    end
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% Reference visibility, C1 or 2, p(correct comparison test stronger) 
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0; 

p.iPlot = p.iPlot+1; 
var = squeeze(mean(val_C_comparison,1)); 
varName = cVarNames{iC}; 
varShortName = 'p(correct reference comparison)';
taskType = 'comparison';
referenceAnnotation = 1; 

group_a.outOfNum = a(1).(varName).OutOfNum';

twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_comparison,referenceAnnotation);

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_referenceComparison',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% RT 
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
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
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
if iC==1 % test stronger
    yticks(0:0.02:0.06)
    ylim([0 0.065])
elseif iC==2 % saw tilt 
    yticks(0:0.05:0.18)
    ylim([0 0.18])
end
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
        ylim([0 0.08])
        yticks(0.:0.02:0.08)
        starText = {'***','***','***','***'}; % main effect, IV, IN, VN
    elseif iC==2
        yPs = [0.88, 0.78, 0.68];
        ylim([0 0.17])
        yticks(0.:0.05:17)
        starText = {'n.s.'}; % '*','**'}; 
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
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-dALLOWPSTRANSPARENCY','-dNOSAFER','-transparent','-p10') 
end

%% AUC of CvA by site 
figure
set(gcf,'Position',[0 0 310 280])

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
ylabel('AUC')
if iC==1 % test stronger
    yticks(0:0.02:0.06)
    ylim([0 0.065])
elseif iC==2 % saw tilt 
    yticks(0:0.05:0.18)
    ylim([0 0.18])
end

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
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-dALLOWPSTRANSPARENCY','-dNOSAFER','-transparent','-p10') 
end

%% SDT: Visibility comparison to reference Detection d'
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0;

p.iPlot = p.iPlot+1;
val_C_dprime = g(iC).comparisonByLL.dprimeDetect.val; % s x ll x att
var = squeeze(mean(val_C_dprime,1)); % average subjects 
varName = cVarNames{iC};
cVarShortNames = {'test stronger'};
varShortName = sprintf('d''\nReference comparison ');
taskType = 'dprime';

referenceAnnotation = 1; 
a_C_dprime.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime,referenceAnnotation);
ylabel({'Comparison \itd''\rm'},'FontSize',p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_C%d(x)_dprime',numel(subjectIDs),s.site,iC);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% SDT: Visibility comparison to reference criterion
figure
set(gcf,'Position',[0 0 380 280])
p.iPlot = 0;

p.iPlot = p.iPlot+1;
val_C_criterion = g(iC).comparisonByLL.criterionDetect.val;
var = squeeze(mean(val_C_criterion,1));
varName = cVarNames{iC};
cVarShortNames = {'test stronger'};
varShortName = sprintf('Criterion\nReference comparison');
taskType = 'criterion';

referenceAnnotation = 1; 
a_C_criterion.(cVarNames{iC}) = twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_criterion,referenceAnnotation);
ylabel({'Comparison criterion'},'FontSize',p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_C%d(x)_criterion',numel(subjectIDs),s.site,iC);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType), '-transparent','-p10')
end

%% Discrimination d'(conditioned on when the stimulus was present)
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

%% Save variables 
% g is n x val 
% g2 is data averaged to group x val 
g2.val_x    = val_x; 
g2.val_logx = val_logx; 
g2.val_C    = val_C; 
g2.val_A    = val_A; 

