function [g] = twcf_plotGroup(PSI,subjectIDs,s,cFieldNames,iC,sites,a)

% 1.4 cue gab dis
% var is structure of group analyses (from behav_analysis) 

%% Setup 
p = twcf_analysisParams;
set(0,'DefaultAxesTitleFontWeight','normal');

% Figure directory 
figDir = sprintf('%s/%s', s.groupFigDir, datestr(now,'yymmdd')); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Get group variables where first index is subject idx 2 is attention 
gFieldNames = {'AUC',...
    'alphaC','betaC','gammaC','lambdaC',...
    'alphaA','betaA','gammaA','lambdaA',...
    'alphaR','betaR','gammaR','lambdaR'}; 

for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
    for iS = 1:numel(PSI)
        for iG = 1:numel(gFieldNames)
            g.(gFieldNames{iG}).val(iS,iAtt) = PSI(iS).(cFieldNames{iC}).(gFieldNames{iG})(iAtt);
        end
        g.pCorrect.val(iS,:,iAtt) = PSI(iS).(cFieldNames{iC}).A.yMean(:,iAtt); 
        g.pSeen.val(iS,:,iAtt) = PSI(iS).(cFieldNames{iC}).(cFieldNames{iC}).yMean(:,iAtt); 
    end
end

%% Plot group vars by att 
for iG = 1:numel(gFieldNames)
    y = g.(gFieldNames{iG}).val;
    x = [1 2 3];
    
    if s.plotFigs
        figure
        figureStyle
        hold on
        set(gcf,'Position',[0 0 200 280])
        title(sprintf('expt: %s\nsite: %s | n = %d', s.exptName, s.site, numel(PSI)))

        for iS = 1:numel(PSI)
            for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
                plot(x,y(iS,:),'Color',p.style.greyColor)
                scatter(x(iAtt),y(iS,iAtt),p.style.szSml,'filled','MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha,'MarkerEdgeColor','w')
            end
        end
        for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
            val = y(:,iAtt); 
            val = val(~isinf(val)); 
            e = errorbar(x(iAtt)+p.style.xOffset,mean(val,'omitnan'),std(val,'omitnan')/sqrt(numel(PSI)),...
                'marker','_','MarkerSize',p.style.sz/5,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerEdgeColor',p.style.attColors(iAtt,:),...
                'CapSize',p.style.errCapSize,'lineWidth',5);
            e.Color = p.style.attColorsMuted(iAtt,:);
        end
        xlim([x(1)-p.style.xOffset*5 x(3)+p.style.xOffset+p.style.xOffset*5])
        ylabel(sprintf('%s',gFieldNames{iG}, cFieldNames{iC}))
        xticks([1,2,3])
        xticklabels({'I','N','V'})
    end

    valid = y(:,3);
    neutral = y(:,2);
    invalid = y(:,1);
    [h pVal ci stats] = ttest(invalid,valid);
    g.(gFieldNames{iG}).h = h; 
    g.(gFieldNames{iG}).pVal = pVal; 
    g.(gFieldNames{iG}).ci = ci; 
    g.(gFieldNames{iG}).stats = stats; 
    g.(gFieldNames{iG}).h = h; 
    
    % Save figures
    if s.saveFigs
        figTitle = sprintf('%s_%s_%s_group_n%d_%s_%s',...
            s.exptFolder,s.site,cFieldNames{iC},...
            numel(PSI),gFieldNames{iG},datestr(now,'yymmdd'));
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    end
end

%% p(seen) and p(correct) by attention 
gFields = {'pCorrect','pSeen'};
for iG = 1:numel(gFields)
    y = g.(gFields{iG}).val; % subjects x contrast x att 
    x = 1:7; 
    
    if s.plotFigs
        figure
        figureStyle
        hold on
        set(gcf,'Position',[0 0 200 280])
        title(sprintf('expt: %s\nsite: %s | n = %d', s.exptName, s.site, numel(PSI)))

        for iS = 1:numel(PSI)
            for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
                plot(x,y(iS,:,iAtt),'Color',p.style.greyColor)
                for iL = 1:7
                    scatter(x(iL),y(iS,iL,iAtt),p.style.szSml,'filled','MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha/5,'MarkerEdgeColor','w')
                end
            end
        end
        for iL = 1:7
            for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
                val = y(:,iL,iAtt);
                val = val(~isinf(val));
                e = errorbar(x(iL)+p.style.xOffset,mean(val,'omitnan'),std(val,'omitnan')/sqrt(numel(PSI)),...
                    'marker','_','MarkerSize',p.style.sz/5,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerEdgeColor',p.style.attColors(iAtt,:),...
                    'CapSize',p.style.errCapSize,'lineWidth',5);
                e.Color = p.style.attColorsMuted(iAtt,:);
            end
        end
        xlim([x(1)-p.style.xOffset*5 x(end)+p.style.xOffset+p.style.xOffset*5])
        xticks(1:7)
        switch gFields{iG}
            case 'pCorrect'
                ylabel('p(correct)')
            case 'pSeen'
                ylabel(sprintf('p(seen) %s', cFieldNames{iC}))
        end
    end

    for iL = 1:7
        valid = y(:,iL,3);
        neutral = y(:,iL,2);
        invalid = y(:,iL,1);
        [h pVal ci stats] = ttest(invalid,valid);
        g.(gFields{iG}).h(iL) = h;
        g.(gFields{iG}).pVal(iL) = pVal;
        g.(gFields{iG}).ci(iL,:) = ci;
        g.(gFields{iG}).stats(iL) = stats;
    end
    
    % Save figures
    if s.saveFigs
        figTitle = sprintf('%s_%s_%s_group_n%d_%s_%s',...
            s.exptFolder,s.site,cFieldNames{iC},...
            numel(PSI),gFields{iG},datestr(now,'yymmdd'));
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    end
end

%% Get group sdt detection vars 
cNames = {'comparison','comparisonByLL'};
for iCN = 1:numel(cNames)
    for iAtt = 1:3 % numel(PSI(1).(cFieldNames{iC}).attConds)
        for iS = 1:numel(PSI) % subjects
            g.(cNames{iCN}).dprimeDetect.val(iS,:,iAtt) = a(iS).sdt.(cNames{iCN}).dprimeDetect(:,iAtt);
            g.(cNames{iCN}).criterionDetect.val(iS,:,iAtt) = a(iS).sdt.(cNames{iCN}).criterionDetect(:,iAtt);
            g.(cNames{iCN}).det_pfa.val(iS,:,iAtt) = a(iS).sdt.(cNames{iCN}).det_pfa(:,iAtt);
            g.(cNames{iCN}).det_ph.val(iS,:,iAtt) = a(iS).sdt.(cNames{iCN}).det_ph(:,iAtt);
            % discrimination sdt
            g.dprimeDiscrim.val(iS,:,iAtt) = a(iS).sdt.dis_correct_V.dprimeDis(:,iAtt);
            g.criterionDiscrim.val(iS,:,iAtt) = a(iS).sdt.dis_correct_V.criterionDis(:,iAtt);
        end
    end
end

%%  Write to table then export csv 
if s.exportcsv
    tableHeaders = {'site', 'att','subjectID',...
        'AUC',...
        'alphaC','betaC','gammaC','lambdaC',...
        'alphaA','betaA','gammaA','lambdaA',...
        'alphaR','betaR','gammaR','lambdaR'};
    count = 1;
    for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
        for iS = 1:numel(PSI)
            V.att(count) = PSI(1).(cFieldNames{iC}).attConds(iAtt);
            V.subjectID(count) = iS;
            switch sites{iS}
                case 'BU'
                    V.site(count) = 0;
                case 'UCI'
                    V.site(count) = 1;
            end
            for iG = 1:numel(gFieldNames)
                V.(gFieldNames{iG})(count) = PSI(iS).(cFieldNames{iC}).(gFieldNames{iG})(iAtt);
            end
            count = count+1;
        end
    end
    % T = struct2table(V);
    sz = [numel(PSI)*numel(PSI(1).(cFieldNames{iC}).attConds) numel(tableHeaders)];

    varTypes = repmat("double",[1 numel(tableHeaders)]);
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

    for iT = 1:numel(tableHeaders)
        T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
    end

    csvName = sprintf('%s_%s_%s_group',...
        s.exptShortName,s.site,cFieldNames{iC}); % saves to csv folder 
    csvPath = sprintf('%s/%s.csv', s.groupCSVDir, csvName);
    writetable(T,csvPath)
end

%% Export pseen pcorrect csv 
clear V
if s.exportcsv
    tableHeaders = {'site', 'att','subjectID',...
        'contrast','pCorrect','pSeen'};
    count = 1;
    for iAtt = 1:numel(PSI(1).(cFieldNames{iC}).attConds)
        for iS = 1:numel(PSI)
            for iL = 1:7
                for iG = 1:numel(gFields)
                    V.contrast(count) = iL;
                    V.(gFields{iG})(count) = g.(gFields{iG}).val(iS,iL,iAtt);

                    V.att(count) = PSI(1).(cFieldNames{iC}).attConds(iAtt);
                    V.subjectID(count) = iS;
                    switch sites{iS}
                        case 'BU'
                            V.site(count) = 0;
                        case 'UCI'
                            V.site(count) = 1;
                    end
                end
                count = count+1;
            end
        end
    end

    % T = struct2table(V);
    sz = [7*numel(PSI)*numel(PSI(1).(cFieldNames{iC}).attConds) numel(tableHeaders)];

    varTypes = repmat("double",[1 numel(tableHeaders)]);
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

    for iT = 1:numel(tableHeaders)
        T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
    end

    csvName = sprintf('%s_%s_%s_group_prob',...
        s.exptShortName,s.site,cFieldNames{iC}); % saves to csv folder 
    csvPath = sprintf('%s/%s.csv', s.groupCSVDir, csvName);
    writetable(T,csvPath)
end

