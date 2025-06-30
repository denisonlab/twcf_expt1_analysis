
% Test RPF 
load('1.3data_subjects.mat'); % load data 

saveFigs = 1; 
figDir = 'figs/HOHO_simRPF'; 
if ~exist(figDir,"dir")
    mkdir(figDir)
end
figType = 'pdf'; 

%% Define Weibull 
paramNames = {'alpha','beta','gamma','lambda'}; 
PF = @PAL_Weibull; 

stimLevelsFine = 0:0.01:0.5; 

% p.colorPurple = [106 90 167]/255; % purple % high awareness l [63 0 125]/255; % high awareness exposure 
% p.colorGreen = [99 188 70]/255; % green % low awareness % [0 128 0]/255; % low awareness exposure 
p.colorPurple = [144 134 178]/255; 
p.colorGreen = [188 211 168]/255; 

%% Make objective - use the group average parameters from the neutral
% condition 

F1(1) = mean( g(1).alphaA.val(:,2), 'omitnan'); 
F1(2) = mean( g(1).betaA.val(:,2), 'omitnan'); 
F1(3) = mean( g(1).gammaA.val(:,2), 'omitnan'); 
F1(4) = mean( g(1).lambdaA.val(:,2), 'omitnan'); 

% Make subjective 
F2(1) = mean( g(1).alphaC.val(:,2), 'omitnan'); 
F2(2) = mean( g(1).betaC.val(:,2), 'omitnan'); 
F2(3) = mean( g(1).gammaC.val(:,2), 'omitnan'); 
F2(4) = mean( g(1).lambdaC.val(:,2), 'omitnan'); 

% Parameter shifts under "high awareness"
F1_prime = F1; 
F2_prime = F2;


% What if we improve performance slightly? 
% F1_prime(1) = 0.15; 

% What if we shift alpha? 
% F2_prime(1) = F2(1)*0.6; 

% What if we increase the slope?
% F2_prime(2) = F2(2)*2; 

% What if we shift up gamma? 
% F2_prime(3) = 0.3; 

% What if we shift down lambda? 
% F2(4) = 0.1; 
% F2_prime(4) = 0; 

% Vertical shift of PF? 

%% Make the weibulls
fit_F1 = PF(F1, stimLevelsFine);
fit_F2 = PF(F2, stimLevelsFine);

fit_F2 = fit_F2 - 0.15; 

fit_F1_prime = PF(F1_prime, stimLevelsFine); 
fit_F2_prime = PF(F2_prime, stimLevelsFine); 

fit_F2_prime = fit_F2+0.1; 
fit_F2_prime(fit_F2_prime>1)=1; 

%% Plot objective 
figure
set(gcf,'Position',[100 100 500 150])

subplot 131
figureStyle 
plot(stimLevelsFine, fit_F1,'LineWidth',1.5,'Color',p.colorGreen)
plot(stimLevelsFine, fit_F1_prime,'LineWidth',1.5,'Color',p.colorPurple)
xlabel('Contrast')
ylabel('p(correct)')
ylim([0.5 1])
yticks(0.5:0.5:1)

subplot 132
figureStyle 
plot(stimLevelsFine, fit_F2,'LineWidth',1.5,'Color',p.colorGreen)
plot(stimLevelsFine, fit_F2_prime,'LineWidth',1.5,'Color',p.colorPurple)
xlabel('Contrast')
ylabel('p(saw stimulus)')
ylim([0 1])
yticks(0:0.5:1)

subplot 133
figureStyle 
plot(fit_F1, fit_F2,'LineWidth',1.5,'Color',p.colorGreen)
plot(fit_F1_prime, fit_F2_prime,'LineWidth',1.5,'Color',p.colorPurple)
xlabel('p(correct)')
ylabel('p(saw stimulus)')
ylim([0 1])
yticks(0:0.5:1)
xlim([0.5 1])
xticks(0.5:0.5:1)

if saveFigs 
    figTitle = sprintf('lambda_shift');
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figType))
end



