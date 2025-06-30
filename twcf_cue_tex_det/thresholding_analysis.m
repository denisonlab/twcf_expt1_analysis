%% load data

saveFigs = 1; 
exptName   = 'cued texture detection';
exptShortName = 'twcf_cue_tex_det'; 
site       = 'BU';
dataFolder = 'thresholding';

% * * * change here * * *
dataFile   = 'twcf_cue_tex_det_BU_S0002_thresholding_20241217_152630';
% currentFolder = '/Users/karen/Documents/github/TWCF_FOHO/twcf_expt1_data_BU'; % lab desk
currentFolder = '/Users/kantian/Dropbox/github/TWCF_FOHO/twcf_expt1_data_BU'; % laptop
% * * * * * *

dataFileSplitName = strsplit(dataFile,'_');
subjectID  = dataFileSplitName{6}; % 'S0006';
date       = dataFileSplitName{8}; % '20220614_153756'; 
time       = dataFileSplitName{9}; % '20220614_153756'; \

subjectFolder = sprintf('%s/%s/%s/%s',currentFolder,exptShortName,subjectID,dataFolder);  
load(sprintf('%s/%s.mat',subjectFolder,dataFile))

addpath (genpath(pwd))

if saveFigs 
    figFolder = sprintf('%s/figs',subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% plot

figure;

sgtitle(sprintf('expt: %s | site: %s \n subject: %s | date: %s', exptName, site, subjectID, und2space(date)))

subplot(1,2,1); hold on;
for i_qID = 1:3
    plot( data.q(i_qID).intensity(1:40) );
end
plot([1,40], p.quest.tGuess*[1,1], 'k--')
xlabel('figure present @ postcue trial #')
ylabel('log_{10} line length (deg)')
title(['prior = ' num2str(p.quest.tGuess) ', final = ' num2str(median(QuestMean(data.q)))])
legend('track 1', 'track 2', 'track 3', 'prior')

subplot(1,2,2); hold on;
for i_qID = 1:3
    plot( 10.^data.q(i_qID).intensity(1:40) );
end
plot([1,40], 10.^p.quest.tGuess*[1,1], 'k--')
xlabel('figure present @ postcue trial #')
ylabel('line length (deg)')
title(['prior = ' num2str(10.^p.quest.tGuess) ', final = ' num2str(median(10.^QuestMean(data.q)))])

if saveFigs 
    figFolder = sprintf('%s/figs',subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
    figTitle = sprintf('thresholding_%s_%s_%s',site,subjectID,date); 
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle)) 
end

%% Get thresholding lengths 
makePlot = 1; 
threshDir = p.setup.threshDir; 
distFromScreen_inCm = 75; 
pixels_perCm = p.screen.pixels_perCm; 
lineLength_inDeg_list = expt_getThresholds(threshDir, distFromScreen_inCm, pixels_perCm, makePlot);