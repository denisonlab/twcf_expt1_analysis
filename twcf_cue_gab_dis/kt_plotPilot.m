function kt_plotPilot

%% Set path
analysisDir = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_analysis_BU';
cd(analysisDir)
addpath(genpath( [analysisDir '/helper_functions'] ))

saveFigs = 1; 
if saveFigs
    figFolder = sprintf('%s/%s/figs',analysisDir,setup.exptName);
    if ~exist(figFolder)
        mkdir(figFolder)
    end
end

%% Load data
% filtered noise 
filename = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_data_BU/twcf_cue_gab_det/ERPilot/piloting/twcf_cue_gab_det_BU_ERPilot_piloting_20230308_090725.mat';

% uniform noise 
filename = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_data_BU/twcf_cue_gab_det/ERPilot2/piloting/twcf_cue_gab_det_BU_ERPilot2_piloting_20230308_141407.mat';

% kt uniform noise, constant noise before targets (0.2) 
filename = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_data_offsite/twcf_cue_gab_det/ktpilot/piloting/twcf_cue_gab_det_offsite_ktpilot_piloting_20230308_163540_trial_98_of_block_2.mat'; 

load(filename)

%% sort data by postcue contrast 
clear correctDis
clear RT
for iContrast = 1:7
    idx = data.contrastID_postcue==iContrast & data.goodTrial==1; 
    correctDis{iContrast,:} = data.correctDis(idx);
    RT{iContrast,:} = data.RT(idx); 
end

%% Plot (56 points per contrast) 
sz = 100; 
color = [0.5 0.5 0.5]; 

figure
subplot 211 
figureStyle 
hold on 
for iContrast = 1:7
    scatter(iContrast,mean(correctDis{iContrast,:}),sz,'MarkerFaceColor',color)
end
ylabel('p(correct discrimination)')
xlim([0 8])
xticks(1:7)
xticklabels(num2cell(p.stim.periph.gabContrast_list))

subplot 212
figureStyle 
hold on 
for iContrast = 1:7
    scatter(iContrast,median(RT{iContrast,:}),sz,'MarkerFaceColor',color)
end
ylabel('Median reaction time (s)')
xlim([0 8])
xticks(1:7)
xticklabels(num2cell(p.stim.periph.gabContrast_list))
xlabel('Grating contrast')

sgtitle(sprintf('%s\n%s | %s | noise type = %s',...
    und2space(setup.exptName),...
    setup.subjectID,...
    setup.site,...
    p.stim.gab.noiseType))

if saveFigs
    figTitle = sprintf('%s_rt_discrim',setup.dataFilename);
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle))
end
