function runBehav
% Behavioral analysis for TWCF_cue_tex_det fmri
% 
% Combines session level-data into subject-level data. Makes basic
% behavioral plots: 1) p(correct), 2) p(seen), 3) p(seen) vs. p(correct),
% 4) AUC of (3).
% Combines subject-level data into group-level.

%% Custom settings

%%% --- Subject IDs ----
% --- BU ----
% main task
subjectIDs_BU = {'S0144'}; % 'S0004'}; % 'S0004','S0097'
% validation block
% subjectIDs_BU = {'S0148','S0152','S0003','S0002','S0150','S0097','S0004','S0147','S0144'}; 

% --- UCI ---
subjectIDs_UCI = {}; 
subjectIDs = [subjectIDs_BU subjectIDs_UCI];

% --- Specify path to TWCF parent folder ---
% s.baseDir       = '/Users/karen/Dropbox/github/TWCF_FOHO'; % lab computer 
s.baseDir       = '/Users/kantian/Dropbox/github/TWCF_FOHO'; % laptop

%% Setup
s.plotFigs      = 1; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 1; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'png'; % 'svg 

s.exptName      = 'cue texture detection'; % 'cued texture detection'
s.exptShortName = 'twcf_cue_tex_det_fmri'; % 'twcf_cue_tex_det' 
s.exptFolder    = 'main_expt'; % 'main_expt' 'validation'

s.site          = 'UCI+BU'; % 'BU' 'UCI'
sites           = [ repmat({'BU'},[1,numel(subjectIDs_BU)]) repmat({'UCI'},[1,numel(subjectIDs_UCI)]) ];  
s.sites         = sites; 

% unlocked = 0; % different directories depending on whether locked or unlocked data repositories are used 

%% Loads and compiles session data to subject level structure 'dataAll'
s.dataRootName  = 'twcf_expt1_data_BU'; 

addpath(genpath(s.baseDir))

s.analDir     = sprintf('%s/twcf_expt1_analysis_public/%s',s.baseDir,s.exptShortName); 
cd(s.analDir)

s.groupFigDir = sprintf('%s/figs',s.analDir);
if ~exist(s.groupFigDir, 'dir')
    mkdir(s.groupFigDir)
end

s.groupCSVDir = sprintf('%s/twcf_expt1_stats_BU/%s/data',s.baseDir,s.exptShortName);
if ~exist(s.groupCSVDir, 'dir')
    mkdir(s.groupCSVDir)
end

%% Compiles data 
clear dataAll
[dataAll,subjectIDs] = twcf_behavCompile_fmri(subjectIDs, s); % sites 
cd(s.analDir)

%% Compile threshold data
s.exptFolder = 'thresholding'; 
[dataAllThresh,subjectIDs] = behav_compileThresh(subjectIDs, s); % sites 
cd(s.analDir)

s.plotFigs = 1;
s.saveFigs = 1; 
s.figType = 'png'; 
plotThreshTiming(dataAllThresh,s)

%% Plot discrim detect by line length, fit curve 
% % * * * change here * * *
% subjectID = 'S0004'; % subject of interest
% % * * * * * * * * * * * *
clear a; clear PSI;
for i = 1:numel(subjectIDs)
    s.dataFolder    = sprintf('%s/twcf_expt1_data_BU',s.baseDir);
    % s.dataFolder    = sprintf('%s/twcf_expt1_data/%s',s.baseDir,sites{i});
    subjectID       = subjectIDs{i};
    [a(i),PSI(i)]   = twcf_behavAnalysis(dataAll,subjectIDs,subjectID,s); % PSI(i)
    % close all % close figs 
end

%% Do group analysis! 
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'svg'; % svg
s.exportcsv     = 0; 

cFieldNames = {'C1','C2','CX'}; % saw figure, saw shape, saw figure not shape
clear g
for iC = 1:numel(cFieldNames)
    g(iC) = twcf_plotGroup(PSI,subjectIDs,s,cFieldNames,iC,sites,a);
end

%% save workspace here %% 
% save('twcf_fmri_data_subjects.mat')

%% Group figs
addpath(genpath(s.baseDir))
cd(s.analDir)
s.plotFigs      = 1; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 1; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'pdf'; % svg
for iC = 1 % :2 % :3 % plot just C as "saw fig" 
    g2 = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g);
end



