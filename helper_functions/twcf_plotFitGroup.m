function a = twcf_plotFitGroup(var,varName,varShortName,a,p,taskType,s,varSubjects,referenceAnnotation,subjectA)
% a = twcf_plotFit(var,varName,a,p,taskType,figureType)
% Plots var (e.g. A(x), C(x) by stim strength (x) by attention condition 
% And fits psychometric function
% Expts 1.1-1.4 
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
if nargin < 9 
    referenceAnnotation = 0; % default no reference annotation 
end
if nargin < 10 
    plotConditionedOnPrecue = 0; 
else 
    plotConditionedOnPrecue = 1; 
end

%% Settings
annotateThreshold = 0; % plots line at threshold 
plotError = 1; % plots error bar ± 1 SEM
plotErrorX = 0; % plots error bar ± 1 SEM stimulus strength
plotShadedError = 0; % error w exporting to pdf? 
plotConditionedOnPrecue = 0; 

%% Palamedes fit parameters setup (move to params file) 
PF = p.fit.PF; % @PAL_Gumbel; % psychometric function type to fit, @PAL_Logistic, @PAL_Gumbel (ie log-Weibull), @PAL_Weibull

% Stimulus feature 
switch s.exptShortName
    case {'twcf_cue_tex_det','twcf_cue_tex_dis'}  
        stimVar = 'uniqueLineLengths_inDeg_postcue';
    case {'twcf_cue_gab_det','twcf_cue_gab_dis'}
        stimVar = 'uniqueContrasts_postcue';
end

if contains(varName,'A') % accuracy
    paramsFree        = p.fit.paramsFreeA; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed
elseif contains(varName,'C') % consciousness
    paramsFree        = p.fit.paramsFreeC;
end

switch taskType
    case {'dprime','criterion','criterion_diff',...
            'dprime_stimulus','dprime_feature',...
            'criterion_stimulus','criterion_feature',...
            'dprime_stimulus_UV','dprime_feature_UV',...
            'criterion_stimulus_UV','criterion_feature_UV',...
            'slope_UV'}
        paramsFree = [1 1 1 1];
end

% Stimulus intensities and search grid
switch p.fit.PFtype
    case 'gumbel'
        stimLevels = log(a.(stimVar)); % log space (log10)

        searchGrid.alpha  = p.fit.searchGrid.alpha; % threshold
        searchGrid.beta   = p.fit.searchGrid.beta; % slope (adjust upper bound)
        searchGrid.lambda = p.fit.searchGrid.lambda; % lapse-rate

    case 'weibull'
        stimLevels = a.(stimVar); % linear space <-- here subtract 1 pix?

        searchGrid.alpha  = p.fit.searchGrid.alpha; % threshold
        searchGrid.beta   = p.fit.searchGrid.beta; % slope (adjust upper bound)
        searchGrid.lambda = p.fit.searchGrid.lambda; % lapse-rate

    case 'nakarushton'
        stimLevels = a.(stimVar);

        searchGrid.alpha  = p.fit.searchGrid.alpha; % threshold
        searchGrid.beta   = p.fit.searchGrid.beta; % slope (adjust upper bound)
        searchGrid.lambda = p.fit.searchGrid.lambda; % lapse-rate
end

stimLevelsFine = min(stimLevels):(max(stimLevels)-min(stimLevels))/1000:max(stimLevels);

%% Plot
% subplot (1,1,1)
hold on

switch taskType
    case 'detection'
        ylims = [0 1]; 
        yticks(0:0.25:1)
        if contains(varName,'A') % accuracy
            ylabel(sprintf('p(%s)', varShortName),'FontSize',p.style.textAxisSize)
        elseif contains(varName,'C') % consciousness
            ylabel(sprintf('p("%s")', varShortName),'FontSize',p.style.textAxisSize)
        end
    case 'discrimination'
        ylims = [0.4 1];
        yticks(0:0.25:1)
        yline(0.5,'--')
        if contains(varName,'A') % accuracy
            ylabel(sprintf('p(%s)', varShortName),'FontSize',p.style.textAxisSize)
        elseif contains(varName,'C') % consciousness
            ylabel(sprintf('p("%s")', varShortName),'FontSize',p.style.textAxisSize)
        end
    case 'comparison'
        ylims = [0 1]; 
        yticks(0:0.25:1)
        ylabel(sprintf('p(correct reference\n comparison)'),'FontSize',p.style.textAxisSize)
    case 'RT'
        ylims = [0 1.5]; % [0 1]; 
        yticks(0:0.5:3)
        ylabel('Reaction time (s)','FontSize',p.style.textAxisSize)
    otherwise 
        ylims = []; 
        ylabel(sprintf('%s', varShortName),'FontSize',p.style.textAxisSize)
end
if ~isempty(ylims)
    ylim(ylims)
end

figureStyle
%% Reference annotation
if referenceAnnotation
    xline(stimLevels(4),'LineWidth',1,'Color',[0.5 0.5 0.5])
    % yl = ylim;
    % xl = xlim;
    % txt = twcf_annotateStats(stimLevels(4),max(yl),'Reference');
end

switch p.fit.PFtype 
    case 'gumbel'
        xlabel('log_{10} (line length) at cued loc (deg)','FontSize',p.style.textAxisSize)
        xlim([stimLevels(1)-p.style.xBuffer stimLevels(end)+p.style.xBuffer])

    case {'weibull','nakarushton'}
        switch s.exptShortName
            case {'twcf_cue_tex_det'}
                set(gca,'Xscale','log')
                xlabel('Line length (°)','FontSize',p.style.textAxisSize)
                xlim([stimLevels(1)-p.style.xBuffer/10 stimLevels(end)+p.style.xBuffer/10])
            case {'twcf_cue_tex_dis'}
                set(gca,'Xscale','log')
                xlabel('Line length (°)','FontSize',p.style.textAxisSize)
                % xlim([stimLevels(1)-p.style.xBuffer/10 stimLevels(end)+p.style.xBuffer/10])
                xlim([0.165 3.25])
            case {'twcf_cue_gab_det'}
                xlabel('Contrast','FontSize',p.style.textAxisSize)
                % xlim([stimLevels(1)-p.style.xBufferContrast stimLevels(end)+p.style.xBufferContrast])
                xlim([-0.01 0.34])
                
                switch taskType
                    case {'dprime','criterion','criterion_diff',...
                            'dprime_stimulus','dprime_feature',...
                            'criterion_stimulus','criterion_feature',...
                            'dprime_stimulus_UV','dprime_feature_UV',...
                            'criterion_stimulus_UV','criterion_feature_UV','slope_UV'}
                        set(gca,'Xscale','log')
                        xlim([0.078 0.34])
                end

            case {'twcf_cue_gab_dis'}
                xlabel('Contrast','FontSize',p.style.textAxisSize)
                xlim([stimLevels(1)-p.style.xBufferContrast stimLevels(end)+p.style.xBufferContrast])
        end
end
xl = get(gca, 'XLim');

switch s.exptShortName
    case {'twcf_cue_tex_det'}
        xticks(stimLevels([1,2,4,6,7]))
        xticklabels(stimLevels([1,2,4,6,7]))
        set(gca,'Xscale','log')
        xOffset = 1.025;  
    case {'twcf_cue_tex_dis'}
        xticks(stimLevels([1,2,4,5,6,7]))
        xticklabels(stimLevels([1,2,4,5,6,7]))
        set(gca,'Xscale','log')
        xOffset = 1.025; 
    case {'twcf_cue_gab_det'}
        xticks(stimLevels([1,2,3,5,7,8]))
        xOffset = 1.015;
    case {'twcf_cue_gab_dis'}
        xticks(stimLevels([1,3,4,5,6,7]))
        xOffset = 1.015;
end
offDiff = xOffset-1;
offsetFixed = [1-offDiff 1 1+offDiff]; 

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
switch taskType
    case 'criterion_diff'
        nAtt = 1; 
    otherwise 
        nAtt = numel(a.attHeaders); 
end
for iAtt = 1:nAtt
    if strcmp(s.exptShortName,'twcf_cue_gab_det') && (strcmp(taskType,'dprime') || ...
            strcmp(taskType,'criterion') || strcmp(taskType,'criterion_diff') || ...
            strcmp(taskType,'dprime_stimulus') || strcmp(taskType,'criterion_stimulus'))
        stimIdx = 2:numel(stimLevels); 
    else 
        stimIdx = 1:numel(stimLevels); 
    end
    for iL = stimIdx
        x = stimLevels(iL);
        y = var(iL,iAtt);
        % y = y(~isnan(y)); 
        % y = mean(y,'omitnan');
        % scatter(x,y,p.style.sz,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',p.style.alpha,'MarkerEdgeColor','w');
        a.(varName).yMean(iL,iAtt) = y; 
        % Format values for Palamedes PF fit
        NumPos(iAtt,iL) = round(y*OutOfNum(iAtt,iL)); % sum(var{i,iAtt},'omitnan');
    end
    switch taskType
        case 'detection'
            searchGrid.gamma = 0:0.1:0.5; % guess-rate
        case 'discrimination'
            searchGrid.gamma = 0.5;
        case {'dprime','criterion','criterion_diff',...
                'dprime_stimulus','dprime_feature',...
                'criterion_stimulus','criterion_feature',...
                'dprime_stimulus_UV','dprime_feature_UV',...
                'criterion_stimulus_UV','criterion_feature_UV','slope_UV'}
            searchGrid.gamma = 0;
        case 'comparison'
            searchGrid.gamma = 0.5; 
    end
    switch taskType
        case 'RT'
            paramsFitted = NaN; 
            LL = NaN; 
            fit = NaN; 
        otherwise 
            % Fit a psychometric function to data using a Maximum Likelihood criterion
            [paramsFitted, LL, exitflag] = PAL_PFML_Fit(stimLevels, NumPos(iAtt,:), OutOfNum(iAtt,:), searchGrid, paramsFree, PF);
            fit = PF(paramsFitted, stimLevelsFine);
    end
    
    switch taskType
        case {'detection','discrimination'}
            
            switch s.exptShortName
                case {'twcf_cue_tex_det','twcf_cue_tex_dis'}
                    stimLevelsFinePlot = log(stimLevelsFine);
                otherwise
                    stimLevelsFinePlot = stimLevelsFine;
            end
            l(iAtt) = plot(stimLevelsFine,fit,'Color',p.style.attColors(iAtt,:),'linewidth',p.style.fitLineWidth);

            if annotateThreshold
                findThreshold = paramsFitted(1) < stimLevelsFine;
                thresholdIdx = find(findThreshold, 1, 'first');

                x = [paramsFitted(1) paramsFitted(1)];   % adjust length and location of arrow
                y = [fit(thresholdIdx) ylims(1)];

                [y_norm] = y_to_norm_v2(y(1),y(2));
                [x_norm] = x_to_norm_v2(x(1),x(2));

                annotation('textarrow',x_norm,y_norm,'FontSize',13,'Linewidth',2,'Color',p.style.attColors(iAtt,:),'LineStyle',':')
            end
    end
    
    switch taskType
        case {'dprime','criterion','comparison','RT'}
            color = p.style.attColors(iAtt,:); 
            plot(stimLevels(stimIdx),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        case 'criterion_diff'
            yline(0,'k')
            color = p.style.attColors(4,:);
            plot(stimLevels(stimIdx),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        % case {'dprime_stimulus'}
        %     color = p.style.attColors(iAtt,:);
        %     plot(stimLevels(stimIdx),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth,'LineStyle','--')
        % case {'dprime_feature'}
        %     color = p.style.attColors(iAtt,:);
        %     plot(stimLevels(stimIdx)*xOffset,var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        % case {'dprime_stimulus_UV'}
        %     color = p.style.attColorsKnit(iAtt,:,1);
        %     plot(stimLevels(stimIdx),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth,'LineStyle','--')
        % case {'dprime_feature_UV'}
        %     color = p.style.attColorsKnit(iAtt,:,1);
        %     plot(stimLevels(stimIdx)*xOffset,var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        % case {'criterion_stimulus_UV'}
        %     color = p.style.attColorsKnit(iAtt,:,1);
        %     plot(stimLevels(stimIdx)*offsetFixed(iAtt),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth,'LineStyle','--')
        % case {'criterion_feature_UV'}
        %     color = p.style.attColorsKnit(iAtt,:,1);
        %     plot(stimLevels(stimIdx)*offsetFixed(iAtt),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        case {'criterion_stimulus','criterion_stimulus_UV','dprime_stimulus','dprime_stimulus_UV'}
            color = p.style.attColorsMutedLight(iAtt,:);
            plot(stimLevels(stimIdx)*offsetFixed(iAtt),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth,'LineStyle',':')
        case {'criterion_feature','criterion_feature_UV','dprime_feature','dprime_feature_UV'}
            color = p.style.attColors(iAtt,:);
            plot(stimLevels(stimIdx)*offsetFixed(iAtt),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        case {'slope_UV'}
            color = p.style.attColors(iAtt,:,1);
            plot(stimLevels(stimIdx)*offsetFixed(iAtt),var(stimIdx,iAtt),'Color',color,'linewidth',p.style.fitLineWidth)
        otherwise
            color = p.style.attColors(iAtt,:);
    end
    
    for iL = stimIdx
        x = stimLevels(iL);
        y = var(iL,iAtt);
        
        switch taskType 
            case {'dprime_stimulus','dprime_feature',...
                    'criterion_stimulus','criterion_feature',...
                    'dprime_stimulus_UV','dprime_feature_UV',...
                    'criterion_stimulus_UV','criterion_feature_UV','slope_UV'}
                % don't plot scatter
            otherwise
                sc = scatter(x,y,p.style.sz,'MarkerFaceColor',color,'MarkerFaceAlpha',1,'MarkerEdgeColor','w','LineWidth',2);
                if strcmp(s.exptShortName,'twcf_cue_gab_det') && iL==1
                    sc.Marker = 'x';
                    sc.MarkerEdgeColor = color;
                    sc.SizeData = p.style.sz*1.4;
                else
                    sc.Marker = 'o';
                end
        end
        if plotError
            condSubjects = varSubjects(:,iL,iAtt);
            se(iL) = std( condSubjects ) / sqrt(numel(condSubjects));
            switch taskType
                case 'criterion_diff'
                    errorbar(x,y,se(iL),'Color',p.style.attColorsMutedLight(4,:),'CapSize',0,'LineWidth',2)
                % case {'dprime_stimulus','dprime_stimulus_UV'}
                %     errorbar(x,y,se(iL),'Color',p.style.attColorsMutedLight(iAtt,:),'CapSize',0,'LineWidth',2)
                % case {'dprime_feature','dprime_feature_UV'} 
                %     errorbar(x*xOffset,y,se(iL),'Color',p.style.attColorsMutedLight(iAtt,:),'CapSize',0,'LineWidth',2)
                case {'criterion_stimulus','criterion_stimulus_UV','dprime_stimulus','dprime_stimulus_UV'}
                    errorbar(x*offsetFixed(iAtt),y,se(iL),'Color',p.style.attColorsMutedLight(iAtt,:),'CapSize',0,'LineWidth',1,'LineStyle','none')
                case {'criterion_feature','criterion_feature_UV','dprime_feature','dprime_feature_UV','slope_UV'}
                    errorbar(x*offsetFixed(iAtt),y,se(iL),'Color',p.style.attColors(iAtt,:),'CapSize',0,'LineWidth',1,'LineStyle','none')
                % case {'slope_UV'}
                %     errorbar(x*offsetFixed(iAtt),y,se(iL),'Color',p.style.attColors(iAtt,:),'CapSize',0,'LineWidth',2)
                otherwise
                    errorbar(x,y,se(iL),'Color',p.style.attColorsMutedLight(iAtt,:),'CapSize',0,'LineWidth',2)
            end
        end

        if plotErrorX
            for iS = 1:numel(condSubjects)
                lls(iS) = a.uniqueLineLengths_inDeg_postcue_S(iS,iL);
            end
            se = std(lls) / sqrt(numel(condSubjects));
            errorbar(x,y,se,'horizontal','Color',p.style.attColorsMutedLight(iAtt,:),'CapSize',0,'LineWidth',2)
        end

    end
    
    % if plotShadedError
    %     clear y
    %     for iL = 1:numel(stimLevels)
    %         y(iL) = var(iL,iAtt) + se(iL); 
    %         yNeg(iL) = var(iL,iAtt) - se(iL); 
    %     end
    %     x = [stimLevels flip(stimLevels)]; 
    %     y = [y flip(yNeg)]; 
    %     % pa = patch(x,y,p.style.attColorsExtraLight(iAtt,:),...
    %     %     'EdgeColor','w','LineWidth',0.1,'FaceAlpha',1); 
    %     % Draw patch 
    %     pa = fill(x,y,p.style.attColorsExtraLight(iAtt,:));
    %     % pa.FaceColor = p.style.attColorsExtraLight(iAtt,:); 
    %     % pa.EdgeColor = [0 0 0]; 
    %     % pa2 = area(stimLevels,yNeg);
    %     % pa2.FaceColor = 'w'; 
    %     % pa2.EdgeColor = [0 0 0]; 
    %     uistack(pa,'bottom')
    %     % uistack(pa2,'bottom')
    % end

    a.(varName).paramsFitted(iAtt,:) = paramsFitted; 
    a.(varName).LL(iAtt,:) = LL; 
    % a.(varName).exitflag(iAtt,:) = exitflag; 
    a.(varName).fit(iAtt,:) = fit; 
end
a.(varName).NumPos = NumPos;
a.(varName).OutOfNum = OutOfNum; 

%% Plot conditioned on precue 
clear dat
fieldNames = {'correctDet_conditionedOnPrecue','correctDis_conditionedOnPrecue'};
stimIdx = 1:numel(stimLevels); 
if plotConditionedOnPrecue
    for iF = 1:numel(fieldNames)
        % Average subjects
        for iS = 1:numel(subjectA)
            for iAtt = 1:numel(a.attHeaders)
                for iL = stimIdx
                    dat(iL,iAtt,iS) = mean(subjectA(iS).(fieldNames{iF}){iL,iAtt},'omitnan');
                end
            end
        end
        a.precue.(fieldNames{iF}) = mean(dat,3,'omitnan'); % average subjects
    end

    % Plot 
    for iAtt = 1:numel(a.attHeaders)
        for iL = stimIdx
            x = stimLevels(iL);
            if contains(varName,'A') % accuracy
                fName        = 'correctDis_conditionedOnPrecue'; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed
            elseif contains(varName,'C') % consciousness
                fName        = 'correctDet_conditionedOnPrecue';
            end
            % Plot discrimation
            y = a.precue.(fName)(iL,iAtt);
            sc = scatter(x,y,p.style.sz,'MarkerFaceColor',p.style.attColors(iAtt,:),'MarkerFaceAlpha',1,'MarkerEdgeColor','w','LineWidth',2);
                sc.Marker = '+';
                sc.MarkerEdgeColor = p.style.attColors(iAtt,:);
                sc.SizeData = p.style.sz*1.4;
        end
    end
end

%% Plot by max distractor contrast 
clear dat yD
fieldNames = {'contrast_precue'} ; %  ,'contrast_distractorsMax'};
stimIdx = 1:numel(stimLevels); 
cGreys = cmocean('-gray',10);
cGreys = cGreys(2:end,:);
if plotConditionedOnPrecue
    for iF = 1:numel(fieldNames)
        % Average subjects
        for iS = 1:numel(subjectA)
            for iAtt = 1:numel(a.attHeaders)
                for iL = stimIdx % target strength
                    for iDL = stimIdx  % distractor strength
                        idxD = subjectA(iS).(fieldNames{iF}){iL,iAtt}==iDL-1;
                        dat{iDL,iL,iAtt,iS,:} = idxD; 
                    end
                end
            end
        end
    end
    
    figure
    % Gather average measure
    for iAtt = 1:numel(a.attHeaders)
        subplot (3,1,iAtt)
        figureStyle
        xlabel('Target contrast')
        for iDL = stimIdx
            for iL = stimIdx
                x = stimLevels(iL);
                xC = stimLevels(iDL);

                if contains(varName,'A') % accuracy
                    fName        = 'correctDis'; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed
                    ylim([0.5 1])
                    ylabel('p(correct)')
                elseif contains(varName,'C') % consciousness
                    fName        = 'sawGrating';
                    ylim([0 1])
                    ylabel('p(seen)')
                end

                for iS = 1:numel(subjectA)
                    y = subjectA(iS).(fName)(iL,iAtt);
                    y = y{1,1};
                    y = y(dat{iDL,iL,iAtt,iS,:}); 
                    y = y(y~=-1);
                    yD(iDL,iL,iAtt,iS) = mean(y,'omitnan');
                end

                % Plot
                yP(iL,iDL) = mean(squeeze(yD(iDL,iL,iAtt,:)),'omitnan'); 
                sc = scatter(x,yP(iL,iDL),p.style.sz*0.8,'MarkerFaceColor',cGreys(iDL,:),'MarkerFaceAlpha',1,'LineWidth',1);
                sc.Marker = 'o';
                sc.MarkerEdgeColor = p.style.attColors(iAtt,:);
                % sc.SizeData = p.style.sz*1.4;
            end
            plot(stimLevels,yP(:,iDL),'Color',cGreys(iDL,:))
        end
    end

end

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

