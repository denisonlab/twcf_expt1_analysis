function runBehav
% Behavioral analysis for TWCF Expt 1.4, cue_gab_dis

% Finds and combines session data to subject level 
% And subjects into group structure 

% BU subject IDs 
subjectIDs_BU = {'S0021',... o1 ++++++
    'S0099',... o2 +++
    'S0094',... o3
    'S0093',... o4
    'S0097',... o5
    'S0068',... o6
    'S0002',... o7
    'S0085',... o8 
    'S0070',... o9
    'S0112',... o10
    'S0108',... o11
    'S0073',... o12
    'S0100',... o13 
    'S0101',... o14 
    'S0004',... o15
    };

% UCI subject IDs
subjectIDs_UCI = {'063',... o1 
    '064',... o2
    '066',... o3
    '068',... o4
    '071',... o5
    '072',... o6
    '073',... o7
    '076',... o8
    '078',... o9
    '079',... o10
    '080',... o11
    '082',... o12
    '083',... o13
    '084',... o14
    '086'}; 

subjectIDs = [subjectIDs_BU subjectIDs_UCI];

%% Setup
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'png'; % 'svg 

s.exptName      = 'cue gabor discrimination'; % 'cued texture detection'
s.exptShortName = 'twcf_cue_gab_dis'; % 'twcf_cue_tex_det' 
s.exptFolder    = 'main_expt'; % 'main_expt'

s.site          = 'UCI+BU'; % 'BU' 'UCI' 'UCI+BU'
sites           = [ repmat({'BU'},[1,numel(subjectIDs_BU)]) repmat({'UCI'},[1,numel(subjectIDs_UCI)]) ]; 
s.sites         = sites; 

%% Loads and compiles session data to subject level structure 'dataAll'
% s.baseDir       = '/Users/karen/Dropbox/github/TWCF_FOHO'; % lab computer 
s.baseDir       = '/Users/kantian/Dropbox/github/TWCF_FOHO'; % laptop
s.dataRootName  = 'twcf_expt1_data_BU'; 
% s.dataRootDir   = sprintf('%s/%s',s.baseDir); 
addpath(genpath(s.baseDir))

s.analDir       = sprintf('%s/twcf_expt1_analysis_BU/%s',s.baseDir,s.exptShortName); 
cd(s.analDir)

s.groupFigDir = sprintf('%s/figs',s.analDir);
if ~exist(s.groupFigDir, 'dir')
    mkdir(s.groupFigDir)
end

s.groupCSVDir = sprintf('%s/twcf_expt1_stats_BU/data',s.baseDir);
if ~exist(s.groupCSVDir, 'dir')
    mkdir(s.groupCSVDir)
end

% Get figure and psychometric function fitting parameters 
p = twcf_analysisParams; 

%% Compiles data 
clear dataAll
[dataAll,subjectIDs] = behav_compile(subjectIDs, s);
cd(s.analDir)

%% Compile threshold data
s.exptFolder = 'thresholding'; 
[dataAllThresh,subjectIDs] = behav_compileThresh(subjectIDs, s); % sites 
cd(s.analDir)

s.plotFigs = 1;
s.saveFigs = 1; 
s.figType = 'png'; 
plotThreshTiming(dataAllThresh,s)

%% Plot discrim detect by contrast, fit curve to subjects
% % * * * change here * * *
% subjectID = 'S0004'; % subject of interest
% % * * * * * * * * * * * *
clear a; clear PSI;
for i = 1:numel(subjectIDs) 
    % if UCI and BU unlocked data, use this data directory 
    if any(strcmp(s.sites,'BU')) && any(strcmp(s.sites,'UCI'))
        s.dataFolder    = sprintf('%s/%s/%s', s.baseDir, s.dataRootName, sites{i});
    else % if BU locked data only 
        s.dataFolder    = sprintf('%s/%s', s.baseDir, s.dataRootName);
    end

    subjectID = subjectIDs{i};

    plotPSI = 1; % no plotting of PSI until data unlocked 
    [a(i),PSI(i)] = twcf_behavAnalysis(dataAll, subjectIDs, subjectID, s, plotPSI); % PSI(i)

    close all
end

%% Do group analysis! 
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'svg'; % svg
s.exportcsv     = 0; 

cFieldNames = {'C1','C2','C3'}; % test stronger, saw tilt (second button press), test stronger | saw tilt 
clear g
for iC = 1:numel(cFieldNames)
    g(iC) = twcf_plotGroup(PSI,subjectIDs,s,cFieldNames,iC,sites,a);
end

%% save workspace here %% 
% save('1.4data_subjects.mat')

% load('1.4data_subjects.mat')

%% Fit to group avg 
addpath(genpath(s.baseDir))
cd(s.analDir)
s.saveFigs      = 1;
s.figType       = 'pdf'; 
for iC = 1 % 2 % 1:2
    [g2] = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g);
end


