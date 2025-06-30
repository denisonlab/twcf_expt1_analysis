
% Fitting SDT model to psychometric function fit to data
% See docx 'TWCF_Psi-SDT_modeling' from Brian M for derivation 

% 0. Load experimental data

% 1. Fit the psychometric function for p(saw stimulus) for given stimulus
% strength (x) and attention condition (A) using MLE

% 2. Fit the psychometric function for p(saw feature) for given stimulus
% strength (x) and attention condition (A) using MLE

% 3. Translate PF to SDT
% See eq (9) and (7) 

%% Settings
saveFigs = 1;
figFolder = 'figs/simMu';
if ~exist(figFolder,'dir')
    mkdir(figFolder)
end

%% Define shared parameters 
clear p_S p_F
clear sigma mu

c_S = 0; % detection criterion of seeing STIMULUS (arbitrary)
c_F = 1; % detection criterion of seeing FEATURE (arbitrary)
        
%% PF -> SDT vars 
for iS = 1:30 % subjects
    for iA = 1:3 % attention
        for iX = 1:8 % stimstrengths numel(a(iS).C1.stimLevels)
        
        %% Define condition
        % iX = 1; % idx of stimulus strength; idx 1 = 0
        % iA = 2; % idx of attention condition; idx 2 = neutral
        % iS = 1; % idx of subject
        
        % Get stimulus strengths (true)
        % iX = 2; 
        x = a(iS).C1.stimLevels(iX);
        
        % Find closest iXFit to the true presented stimulus strengths
        iXFit = numel(a(iS).C1.stimLevelsFine(a(iS).C1.stimLevelsFine<x))+1; 
        % disp(iXFit)
        xFit = a(iS).C1.stimLevelsFine(iXFit);
        
        %% Calculate p(seen)
        p_S(iS,iX,iA) = a(iS).C1.fit(iA,iXFit); % Step 1. PF fit of seeing STIMULUS at a given stimulus strength and attention
        p_F(iS,iX,iA) = a(iS).C2.fit(iA,iXFit); % Step 2. PF fit of seeing FEATURE at a given stimulus strength and attention
        
        %% Equation (9) - Sigma
        sigma(iS,iX,iA) = c_S - c_F / ( norminv(1-p_S(iS,iX,iA)) - norminv(1-p_F(iS,iX,iA)) );
        
        %% Equation (7) - Mu
        mu(iS,iX,iA) = c_S - sigma(iS,iX,iA) * norminv(1-p_S(iS,iX,iA));
        
        end
    end
end

%% Instead of fitting to the data -- let us observe the systematic impact of
% changes in mu and sigma, which we expect attention to affect, 
% on performance and visibility 
for iS = 1:30 % subjects
    for iA = 1:7 % attention
        for iX = 1:numel(a(iS).C1.stimLevels)
            sigma_simAtt(iS,iX,iA) = sigma(iS,iX,1);
            mu_simAtt(iS,iX,iA) = mu(iS,iX,1)*iA; % for now, no changes to mu
        end
    end
end

sigma = sigma_simAtt; 
mu = mu_simAtt; 

%% Plot group average SDT distributions and PF 
stimStrengths = squeeze(median(mu,1)); % median bc of outliers? mean 
sigmas = squeeze(median(sigma,1)); 
nAtt = size(mu,3);

% p = simulation_params(stimStrengths); % Get fit parameters 

% sigma_noise = 1; % assume equal variances
% sigma_signal = ones([1 numel(stimStrengths)]); % for now let us assume equal variances across stimulus strengths
% 
% mu_noise = 0; % detection task 
% mu_signal = ones([1 numel(stimStrengths)]).*stimStrengths; 
% 
% c = 0.5; % absolute criterion 
% c_sigma = 0.2; % criterion noise 

% an even prettier color map 
cmap = cmocean('ice',size(stimStrengths,1)+1); 
% cmap = cmap(2:end,:); 
cmap = flip(cmap,1); 

% colors for criterion
% pale red, classic red 
cRed = [255 127 127;...
    255 0 0]/255; 

% pd_noise = makedist('Normal','mu',mu_noise,'sigma',sigma_noise); 
for iX = 1:size(stimStrengths,1)
    for iA = 1:3
        pd_signal(iX,iA) = makedist('Normal','mu',stimStrengths(iX,iA),'sigma',sigmas(iX,iA)); 
    end
end

%% Plot SDT distributions (1D)
figure
% set(gcf,'Position',[100 100 700 300])
% pl(1) = plot(pd_noise); 
% pl(1).Color = [0.5 0.5 0.5]; 
% pl(1).LineStyle = '--'; 
% pl(1).LineWidth = 1.5;
for iA = 1:3
    subplot (3,1,iA)
    figureStyle
    hold on
    for iX = 1:size(stimStrengths,1)
        pd = pd_signal(iX,iA);
        pl(iX) = plot(pd);
        pl(iX).Color = cmap(iX+1,:);
        if iX==1
            pl(iX).LineStyle = '--'; 
            pl(iX).Color = [0.5 0.5 0.5]; 
        end
        if iA==1
            % pl(iX).LineStyle = ':';
            title('Invalid')
            xlabel('')
        elseif iA==2
            % pl(iX)fig.LineStyle = '--';
            title('Neutral')
            xlabel('')
        elseif iA==3
            title('Valid')
            xlabel('Internal signal strength')
        end
        xlim([-4 4])
        ylim([0 1])
        ylabel('pdf')

    end
    xline(c_S,'Color',cRed(1,:),'LineWidth',1.5,'label',{'c_S'},'LabelOrientation','horizontal');
    xline(c_F,'Color',cRed(2,:),'LineWidth',1.5,'label',{'c_F'},'LabelOrientation','horizontal');
end

if saveFigs
    figTitle = sprintf('%s_SDTdistributions_1D',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot SDT distributions (2D)
figure
hold on 
figureStyle
axis square 
axis equal

for iA = 2 % neutral
    
    for iX = 1:size(stimStrengths,1)
        % Draw 2D distributions (along x)
        xC = pd_signal(1,iA).mu;
        yC = pd_signal(iX,iA).mu;
        sigma = pd_signal(iX,iA).sigma; 
        if iX==1
            viscircles([xC,yC],sigma,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1)
        else
        viscircles([xC,yC],sigma,'Color',cmap(iX+1,:),'LineWidth',1)
        end

        % Draw 2D distributions (along y)
        xC = pd_signal(iX,iA).mu;
        yC = pd_signal(1,iA).mu;
        sigma = pd_signal(iX,iA).sigma; 
        if iX==1
            viscircles([xC,yC],sigma,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1)
        else
        viscircles([xC,yC],sigma,'Color',cmap(iX+1,:),'LineWidth',1)
        end

    end
end

% Reference lines 
% --- Noise mean 
xline(pd_signal(1,iA).mu,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
yline(pd_signal(1,iA).mu,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
% --- Criteria (detection)
xline(c_S,'Color',cRed(1,:),'LineWidth',1.5,'label',{'c_S'},'LabelOrientation','horizontal');
xline(c_F,'Color',cRed(2,:),'LineWidth',1.5,'label',{'c_F'},'LabelOrientation','horizontal');
yline(c_S,'Color',cRed(1,:),'LineWidth',1.5,'label',{'c_S'},'LabelOrientation','horizontal');
yline(c_F,'Color',cRed(2,:),'LineWidth',1.5,'label',{'c_F'},'LabelOrientation','horizontal');
% --- Criteria (detection)
c_d = refline(1,0); 
c_d.Color = 'k'; 
text(2.8,2.6,{'c_d'})

% Axes 
xlim([-2.5 3])
ylim([-2.5 3])
xticks(-4:1:4)
yticks(-4:1:4)
xlabel('Internal signal strength (+45°)')
ylabel('Internal signal strength (-45°)')

if saveFigs
    figTitle = sprintf('%s_SDTdistributions_2D',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot mu and sigma by stim strength and att
measures = {'mu','sigma'};

logOn = 1; 

% Real data 
if logOn
    data_x = log(mean(g2.val_x,1)); 
else
    data_x = mean(g2.val_x,1); 
end

figure
set(gcf,'Position',[100 100 600 200])
for iM = 1:2
    subplot (1,2,iM)
    figureStyle
    hold on
    
    % if iA==1
    %     title('Invalid')
    % elseif iA==2
    %     title('Neutral')
    % else
    %     title('Valid')
    % end
    if logOn
        xlabel('log(stimulus strength)')
    else
        xlabel('Stimulus strength')
    end
    ylabel(measures{iM})

    for iA = 1:3
        if iM==1
            % mu
            plot(data_x,stimStrengths(:,iA)','Color',p.style.attColors(iA,:))
        elseif iM==2
            % sigma
            plot(data_x,sigmas(:,iA)','Color',p.style.attColors(iA,:))
        end
    end
    
    for iA = 1:3
        for iX = 1:8
            scatter(data_x(iX),pd_signal(iX,iA).(measures{iM}), p.style.sz*0.6, ...
                'filled','MarkerFaceColor',p.style.attColors(iA,:),'LineWidth',1,...
                'MarkerEdgeColor','w');
        end
    end

end

if saveFigs
    if logOn 
        figTitle = sprintf('%s_muSigma_Att_SDT_log',a(1).exptShortName);
    else
        figTitle = sprintf('%s_muSigma_Att_SDT',a(1).exptShortName);
    end
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Generate the PF again via simulation from the PDs 
lambda = 0; % lapse rate 
beta = 0; % 0.5; % guess rate on lapse trials 
nTrials = 1000;

for iA = 1:3
    for iX = 1:size(stimStrengths,1)
        noise.mu = stimStrengths(1,iA); % assume neutral mu for now
        noise.sigma = sigmas(1,iA);

        signal.mu = stimStrengths(iX,iA);
        signal.sigma = sigmas(iX,iA);
        
        % saw stimulus
        sdtSim(iA).pS(iX) = simulateObserver(noise, signal, c_S, lambda, beta, nTrials);

        % saw feature
        sdtSim(iA).pF(iX) = simulateObserver(noise, signal, c_F, lambda, beta, nTrials);
    end
end

%% Plot the simulated PF
figure
hold on 
figureStyle
set(gcf,'Position',[100 100 500 350])

% Real data 
data_x = mean(g2.val_x,1); 

clear pS pF
pS = NaN(3,8); 
pF = NaN(3,8); 
var = 'pSeen'; 
for iA = 1:3
    for iX = 1:8
        pS(iA,iX) = sdtSim(iA).pS(iX).SDTvars.hits/500; % pSeen 
        pF(iA,iX) = sdtSim(iA).pF(iX).SDTvars.hits/500;
    end
end

for iA = 1:3  
    plot(data_x,pS(iA,:),'Color',p.style.attColors(iA,:),'LineWidth',1,'LineStyle','--')
    plot(data_x,pF(iA,:),'Color',p.style.attColors(iA,:),'LineWidth',1.5)
end

% Simulated PF from SDT PD
for iA = 1:3
    for iX = 1:size(stimStrengths,1)
        scatter( data_x(iX), sdtSim(iA).pF(iX).pSeen, p.style.szSml, ...
            'filled','MarkerFaceColor',p.style.attColors(iA,:),...
            'MarkerEdgeColor','w');

        scatter( data_x(iX), sdtSim(iA).pS(iX).pSeen, p.style.szSml, ...
            'filled','MarkerFaceColor',p.style.attColors(iA,:),...
            'MarkerEdgeColor','w');
    end
end

xlim([data_x(1)-p.style.xBuffer data_x(end)+p.style.xBuffer])
xticks(data_x)
xtickformat('%0.2f')
xtickangle(30)

xlabel('Contrast')
ylabel('p(seen)')

% Legend
% ax = gca; 
alignment = 'left'; 
kt_annotateStats(0,0.95,'- - saw stimulus',alignment); 
kt_annotateStats(0,0.90,'— saw feature',alignment); 

if saveFigs
    figTitle = sprintf('%s_SDTdistributions_simulatePF',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Simulate objective performance
% The pattern across x and att should match dis d' bc the orientation 
% decision is modeled as bias free

%% 
% p.fit.PFtype = 'gumbel'; % NakaRushton 
% 
% % Fit a psychometric function
% switch p.fit.PFtype
%     case {'gumbel','NakaRushton'}
%         stimLevels = log(data_x); % log space (log10)
%     case 'weibull'
%         stimLevels = data_x; % linear space
% end
% stimLevelsFine = min(stimLevels):(max(stimLevels)-min(stimLevels))/1000:max(stimLevels);
% 
% % Define search grid
% switch p.fit.PFtype
%     case {'gumbel','weibull'}
%         if contains(varName,'A') % accuracy
%             paramsFree        = p.fit.paramsFreeA; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed
%         elseif contains(varName,'C') % consciousness
%             paramsFree        = p.fit.paramsFreeC;
%         end
% 
%         searchGrid.alpha  = p.fit.searchGrid.alpha; % threshold
%         searchGrid.beta   = p.fit.searchGrid.beta; % slope (adjust upper bound)
%         searchGrid.lambda = p.fit.searchGrid.lambda; % lapse-rate
%     case 'NakaRushton'
%         if contains(varName,'A') % accuracy
%             paramsFree        = [1 1 ]; % [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed
%         elseif contains(varName,'C') % consciousness
%             paramsFree        = p.fit.paramsFreeC;
%         end
% 
%         searchGrid.c50 = 1; % should be between their lowest and highest x
%         searchGrid.Rmax = 50:100; % is this in percent or proportion? does it matter? 
%         searchGrid.n = 1; 
%         searchGrid;J = 1; 
% end
% 
% c50 = []; Rmax = []; n = []; M = []; 
% 
% NumPos = a(1).C1.NumPos; 
% OutOfNum = a(1).C1.OutOfNum; 
% 
% idx = 2:8; 

%% 
% figure
% figureStyle
% hold on
% for iAtt = 1
%     p.fit.PF = '@PAL_Gumbel';
%     [paramsFitted, LL, exitflag] = PAL_PFML_Fit(stimLevels(idx), NumPos(iAtt,idx), OutOfNum(iAtt,idx), searchGrid, paramsFree, p.fit.PF);
%     PF = p.fit.PF;
% 
%     fit = PF(paramsFitted, stimLevelsFine);
%     l(iAtt) = plot(stimLevelsFine,fit,'Color',p.style.attColors(iAtt,:),'linewidth',1);
%     xline(paramsFitted(1),'--','Color',p.style.attColors(iAtt,:),'linewidth',1.5) % threshold
% end

%% Calculate detection and discrimination sensitivity 
for iA = 1:3
    for iX = 1:size(stimStrengths,1)
        noise.mu = stimStrengths(1,iA); % assume neutral mu for now
        noise.sigma = sigmas(1,iA);

        signal.mu = stimStrengths(iX,iA);
        signal.sigma = sigmas(iX,iA);
        
        % detection sensitivity
        d_a(iX,iA) = (signal.mu-noise.mu)/sqrt((noise.sigma^2+signal.sigma^2)/2);

        % discrimination sensitivity
        d_d(iX,iA) = (signal.mu*sqrt(2))/signal.sigma;
    end
end

%% Plot detection sensitivity 
logOn = 1; 

% Real data 
if logOn
    data_x = log(mean(g2.val_x,1)); 
else
    data_x = mean(g2.val_x,1); 
end

figure
set(gcf,'Position',[100 100 300 200])
    figureStyle
    hold on

    if logOn
        xlabel('log(stimulus strength)')
    else
        xlabel('Stimulus strength')
    end
    ylabel('Detection sensitivity')

    for iA = 1:3
        plot(data_x,d_a(:,iA)','Color',p.style.attColors(iA,:))
    end
    
    for iA = 1:3
        for iX = 1:8
            scatter(data_x(iX),d_a(iX,iA), p.style.sz*0.6, ...
                'filled','MarkerFaceColor',p.style.attColors(iA,:),'LineWidth',1,...
                'MarkerEdgeColor','w');
        end
    end

if saveFigs
    if logOn 
        figTitle = sprintf('%s_SDT_d_a_log',a(1).exptShortName);
    else
        figTitle = sprintf('%s_SDT_d_a',a(1).exptShortName);
    end
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot discrimination sensitivity 
logOn = 1; 

% Real data 
if logOn
    data_x = log(mean(g2.val_x,1)); 
else
    data_x = mean(g2.val_x,1); 
end

figure
set(gcf,'Position',[100 100 300 200])
    figureStyle
    hold on

    if logOn
        xlabel('log(stimulus strength)')
    else
        xlabel('Stimulus strength')
    end
    ylabel('Discrimination sensitivity')

    for iA = 1:3
        plot(data_x,d_d(:,iA)','Color',p.style.attColors(iA,:),...
            'Marker','o','MarkerSize',8)
    end
    % '-s','MarkerSize',10,...
    % 'MarkerEdgeColor','red',...
    % 'MarkerFaceColor',[1 .6 .6]

    % for iA = 1:3
    %     for iX = 1:8
    %         scatter(data_x(iX),d_d(iX,iA), p.style.sz*0.6, ...
    %             'filled','MarkerFaceColor',p.style.attColors(iA,:),'LineWidth',1,...
    %             'MarkerEdgeColor','w');
    %     end
    % end

if saveFigs
    if logOn 
        figTitle = sprintf('%s_SDT_d_d_log',a(1).exptShortName);
    else
        figTitle = sprintf('%s_SDT_d_d',a(1).exptShortName);
    end
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot detection sensitivity as function of discrimination sensitivity
figure
set(gcf,'Position',[100 100 300 250])
figureStyle
hold on

xlabel('Discrimination sensitivity')
ylabel('Detection sensitivity')

for iA = 1:3
    plot(d_d(:,iA)',d_a(:,iA)','Color',p.style.attColors(iA,:),...
        'Marker','o','MarkerSize',8)
end

% for iA = 1:3
%     for iX = 1:8
%         scatter(d_d(iX,iA),d_a(iX,iA), p.style.sz*0.6, ...
%             'filled','MarkerFaceColor',p.style.attColors(iA,:),'LineWidth',1,...
%             'MarkerEdgeColor','w');
%     end
% end

if saveFigs
    figTitle = sprintf('%s_SDT_detection_vs_discrimination',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot p(seen) as a function of discrimination sensitivity [p(correct)]
figure
set(gcf,'Position',[100 100 300 250])
figureStyle
hold on

xlabel('Discrimination sensitivity')
ylabel('p(saw stimulus)')

for iA = 1:3
    pl = plot(d_d(:,iA)',pS(iA,:),'Marker','o','Color',p.style.attColors(iA,:),'LineWidth',1);
    pl.MarkerSize = 10; 
    pl.MarkerFaceColor = p.style.attColors(iA,:); 
    pl.MarkerEdgeColor = 'w'; 
end

if saveFigs
    figTitle = sprintf('%s_SDT_pSawStimulus_vs_discrimination',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end

%% Plot p(saw feature) as a function of discrimination sensitivity [p(correct)]
figure
set(gcf,'Position',[100 100 300 250])
figureStyle
hold on

xlabel('Discrimination sensitivity')
ylabel('p(saw feature)')

for iA = 1:3
    pl = plot(d_d(:,iA)',pF(iA,:),'Marker','o','Color',p.style.attColors(iA,:),'LineWidth',1);
    pl.MarkerSize = 10; 
    pl.MarkerFaceColor = p.style.attColors(iA,:); 
    pl.MarkerEdgeColor = 'w'; 
end

if saveFigs
    figTitle = sprintf('%s_SDT_pSawFeature_vs_discrimination',a(1).exptShortName);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end








