clear

%% load data

saveFigs = 1; 
exptName   = 'cued texture discrimination';
exptShortName = 'twcf_cue_tex_dis'; 
site       = 'BU';
dataFolder = 'thresholding';

% * * * change here * * *
dataFile   = 'twcf_cue_tex_dis_BU_S0003_thresholding1_20220722_234717';
currentFolder = '/Users/karen/Documents/GitHub/twcf_expt1_data_BU'; % lab desk
% currentFolder = '/Users/kantian/Dropbox/gitHub/twcf_expt1_data_BU'; %
% laptop
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
    l = plot( data.q_cll(i_qID).intensity(1:40) );
    l.Color(4) = 0.75; 
end
plot([1,40], p.quest_cll.tGuess*[1,1], 'k--')
xlabel('trial #')
ylabel('log_{10} center line length (deg)')
title(['prior = ' num2str(p.quest_cll.tGuess) ', final = ' num2str(median(QuestMean(data.q_cll)))])
legend('track 1', 'track 2', 'track 3', 'prior')

subplot(1,2,2); hold on;
for i_qID = 1:3
    l = plot( 10.^data.q_cll(i_qID).intensity(1:40) );
    l.Color(4) = 0.75; 
end
plot([1,40], 10.^p.quest_cll.tGuess*[1,1], 'k--')
xlabel('trial #')
ylabel('center line length (deg)')
title(['prior = ' num2str(10.^p.quest_cll.tGuess) ', final = ' num2str(median(10.^QuestMean(data.q_cll)))])

if saveFigs 
    figFolder = sprintf('%s/figs',subjectFolder);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
    figTitle = sprintf('thresholding_ecc_%s_%s_%s',site,subjectID,date); 
    saveas(gcf,sprintf('%s/%s.png', figFolder, figTitle)) 
end
