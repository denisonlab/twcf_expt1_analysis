
% Load 1.3 grating detection data 
% Load 1.1 texture data

%% Compile sdt measures into sdt structure
% Settings
figDir = 'figs';
figType = 'pdf';
exptName = 'det_gab'; % 'det_gab' 'det_tex'

%% Basic setup 
switch exptName
    case 'det_tex'
        visMeasures = {'sawFigure','sawShape'};
        cIdx = 1:7;
        load('/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_analysis_BU/twcf_cue_tex_det/1.1data_subjects.mat')
    case 'det_gab'
        visMeasures = {'sawGrating','sawOri'};
        cIdx = 2:8; 
        load('/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_analysis_BU/twcf_cue_gab_det/1.3data_subjects.mat')
    otherwise 
        error('expt name unrecognized')
end
p = twcf_analysisParams; 

%% Pull out the relevant vars
clear nh nsignal nnoise nfa h fa z_h z_fa
for iS = 1:30
    % clear sdt
    sdt = a(iS).sdt;
    for iM = 1:2 % saw stimulus, saw grating
        for iAtt = 1:3
            for iC = cIdx
                switch exptName
                    case 'det_gab'
                        iC2 = iC-1; 
                    case 'det_tex'
                        iC2 = iC; 
                end

                loglinear = 1; 

                if loglinear
                    nh(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nh(iC,iAtt) + 0.5; 
                    nfa(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nfa(iC,iAtt) + 0.5;
                    nsignal(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nsignal(iC,iAtt) + 1;
                    nnoise(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nnoise(iC,iAtt) + 1;
                else
                    nh(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nh(iC,iAtt);
                    nfa(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nfa(iC,iAtt);
                    nsignal(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nsignal(iC,iAtt);
                    nnoise(iS,iC,iAtt,iM) = sdt.(visMeasures{iM}).det_nnoise(iC,iAtt);
                end

                % calculate rates
                h(iS,iC,iAtt,iM) = nh(iS,iC,iAtt,iM)/nsignal(iS,iC,iAtt,iM); 
                fa(iS,iC,iAtt,iM) = nfa(iS,iC,iAtt,iM)/nnoise(iS,iC,iAtt,iM); 
                
                % z-score rates 
                z_h(iS,iC,iAtt,iM) = norminv(h(iS,iC,iAtt,iM) );
                z_fa(iS,iC,iAtt,iM) = norminv(fa(iS,iC,iAtt,iM) );
                               
            end
        end
    end
end

%% Check the rates 
h_mean = mean(h(:,:,2,1)); 
h_mean = mean(h(:,:,2,2)); 

fa_mean = mean(fa(:,:,2,1)); 
fa_mean = mean(fa(:,:,2,2));

%% Calculate UV SDT ROC
% plot normalized hit vs. far
% by attention and stimulus strength (single subjects) 
figure
figureStyle
iS = 20; % change subject here
for iC = cIdx
    for iAtt = 1:3
        plot([z_fa(iS,iC,iAtt,1), z_fa(iS,iC,iAtt,2)], [z_h(iS,iC,iAtt,1),z_h(iS,iC,iAtt,2)],...
            'Color',p.style.attColors(iAtt,:),'LineWidth',1.5)
    end
end
% xlim([-2.5 2.5])
% ylim([-2.5 2.5])
refline(1,0)

ylabel('z(h)')
xlabel('z(fa)')
axis equal

%% We have a graded subjective measure with two responses (saw stimulus & saw orientation)
% which we can leverage to fit a line
% to estimate s for the UV SDT model

clear gsdt dprime criterion d_a c_a
switch exptName
    case 'det_gab'
        csIdx = 1:8;
    case 'det_tex'
        csIdx = 1:7;
end
for iS = 1:30
    for iC = csIdx
        for iAtt = 1:3
            
            syms c b % define symbolic
            eqn_s = z_h(iS,iC,iAtt,1) == b * z_fa(iS,iC,iAtt,1) + c; % saw stimulus
            eqn_f = z_h(iS,iC,iAtt,2) == b * z_fa(iS,iC,iAtt,2) + c; % saw feature

            sympref('FloatingPointOutput',true);
            S = solve([eqn_s; eqn_f],[c,b]);

            mu_s = S.c/S.b;
            bV = double(S.b);
            switch exptName
                case 'det_gab'
                    if bV==0 || iC==1
                        s = 0;
                    else
                        s = 1/bV;
                    end
                case 'det_tex'
                    if isempty(bV)
                        s = 1;
                    elseif bV==0
                        s = 1; 
                    else
                        s = 1/bV;
                    end
            end

            % save
            gsdt.s(iS,iC,iAtt) = double(s); 

            for iM=1:2
                % loglinear = 1;
                % [dprime, criterion, d_a, c_a] = kt_dprime(nh,nfa,nsignal,nnoise,loglinear,s);

                % Vanilla SDT
                dprime = z_h(iS,iC,iAtt,iM) - z_fa(iS,iC,iAtt,iM);
                criterion = -0.5*(z_h(iS,iC,iAtt,iM)+z_fa(iS,iC,iAtt,iM));

                % Unequal variance SDT
                d_a = sqrt(2/(1+s^2)) * (z_h(iS,iC,iAtt,iM) - s*z_fa(iS,iC,iAtt,iM));
                c_a = - (sqrt(2)*s / (sqrt(1+s^2)*(1+s)) ) * (z_h(iS,iC,iAtt,iM) + z_fa(iS,iC,iAtt,iM));

                % save
                gsdt.d_a(iS,iC,iAtt,iM) = double(d_a);
                gsdt.c_a(iS,iC,iAtt,iM) = double(c_a);
                gsdt.dprime(iS,iC,iAtt,iM) = dprime;
                gsdt.criterion(iS,iC,iAtt,iM) = criterion;
            end

        end
    end
end

% save sdt 
% save('1.1sdt.mat','gsdt')

%% zero pad (only for gratings) 
fields = fieldnames(gsdt); 
for iF = 2:numel(fields)
    for iS  = 1:30
        for iAtt = 1:3
            for iM = 1:2  
                val = gsdt.(fields{iF})(iS,:,iAtt,iM); 
                gsdt2.(fields{iF})(iS,:,iAtt,iM) = [NaN val]; 
            end
        end
    end 
end

%% Plot to check the differences 
% average subjects
iM = 1; % stimulus vis only 
switch exptName
    case 'det_gab'
        pX = 0:7; 
    case 'det_tex'
        pX = 1:7;
end
figure
% dprime 
subplot 121
figureStyle
for iAtt = 1:3
    x = pX; 
    % UV 
    y = gsdt.d_a(:,:,iAtt,iM); 
    plot(x,mean(y,1,'omitnan'),'-o','Color',p.style.attColors(iAtt,:))
    % Vanilla
    y = gsdt.dprime(:,:,iAtt,iM); 
    plot(x,mean(y,1,'omitnan'),'--x','Color',p.style.attColors(iAtt,:))
end
xlabel('Stimulus strength level')
ylabel('d_{a}')

% criterion
subplot 122
figureStyle
for iAtt = 1:3
    % UV 
    y = gsdt.c_a(:,:,iAtt,iM); 
    plot(x,mean(y,1,'omitnan'),'-o','Color',p.style.attColors(iAtt,:))
    % Vanilla
    y = gsdt.criterion(:,:,iAtt,iM); 
    plot(x,mean(y,1,'omitnan'),'--x','Color',p.style.attColors(iAtt,:))
end
xlabel('Stimulus strength level')
ylabel('c_{a}')

%% 
% syms a b % define symbolic 
% eqn_s = norminv(h) == b * norminv(fa) + a; % saw stimulus 
% eqn_f = norminv(h) == b * norminv(fa) + a; % saw feature 
% 
% mu_s = a/b; 
% 
% s = 1/b; 

%% Calculate sdt measures 
% From Li, Lau & Odegaard 2018 
% d_a = sqrt(2)*mu_s / (sqrt(1+s^2)); 
% c_a = - (sqrt(2)*s / (sqrt(1+s^2)*(1+s)) ) * (norminv(h) + norminv(fa)); 

%% Plot gsdt structure
% Detection d' (overlay stimulus and feature detection) 
% a = twcf_plotFitGroup(var,varName,varShortName,a,p,taskType,s,varSubjects,referenceAnnotation,subjectA); 
iC = 1; 
cVarNames = {'C1','C2'}; 
group_a.attConds = a(1).attConds; 
group_a.attHeaders = a(1).attHeaders; 
group_a.uniqueContrasts_postcue = mean(val_x,1); 
group_a.log10_uniqueContrasts_postcue = mean(val_logx,1); 

% Stimulus detection
val_C_dprime_stimulus = gsdt.dprime(:,:,:,1); % g(1).dprimeDetect.val; % 30 x 8 x 3 
val_C_criterion_stimulus = gsdt.criterion(:,:,:,1);
% Feature detection
val_C_dprime_feature = gsdt.dprime(:,:,:,2);
val_C_criterion_feature = gsdt.criterion(:,:,:,2);

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
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime_stimulus);

%% plot feature
var = squeeze(mean(val_C_dprime_feature,1));
varName = cVarNames{iC};
cVarShortNames = {'saw grating','saw orientation'};
varShortName = sprintf('Detection d''\n%s',cVarShortNames{iC});
taskType = 'dprime_feature';
twcf_plotFitGroup(var,varName,varShortName,group_a,p,taskType,s,val_C_dprime);
ylabel({'Detection \itd''\rm'},'FontSize',p.style.textAxisSize)

if s.saveFigs
    figTitle = sprintf('group_n%d_%s_dprime_overlaid',numel(subjectIDs),s.site);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType),'-transparent','-p10')
end

%% What were the s by attention? 
figure
set(gcf,'Position',[100 100 300 500])
for iAtt = 1:3
    subplot (3,1,iAtt)
    figureStyle
    xticks(cIdx)
    xlim([cIdx(1)-1 cIdx(end)+1])
    for iC = cIdx
        val = gsdt.s(:,iC,iAtt); 
        val = val(val~=0);
        valMean = mean(val,'omitnan'); 
        valSD = std(val,'omitnan')/sqrt(30); 
        errorbar(iC,valMean,valSD,'k')
        bP = bar(iC,valMean); 
        bP.BaseValue = 1; 
        bP.FaceColor = p.style.attColorsMuted(iAtt,:);
        % scatter(iC,val)
    end
end
ylabel('UV SDT s')
xlabel('Contrast level')

if saveFigs
    % figType = 'pdf'; 
    figTitle = sprintf('group_cue_%s_UVsdt_s',exptName);
    export_fig(gcf,sprintf('%s/%s.%s', figDir, figTitle, figType),'-transparent','-p10')
end