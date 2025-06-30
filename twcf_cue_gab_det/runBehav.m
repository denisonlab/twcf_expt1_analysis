function runBehav
% Behavioral analysis for TWCF Expt 1.3, cue_gab_det

% Finds and combines session data to subject level 
% And subjects into group structure 

% BU subject IDs 
subjectIDs_BU = {'S0021',... o1 
    'S0049',... o2 
    'S0093',... o3
    'S0094',... o4 
    'S0009',... o5
    'S0067',... o6
    'S0082',... o7
    'S0068',... o8
    'S0066',... o9
    'S0063',... o10
    'S0091',... o11
    'S0087',... o12
    'S0092',... o13
    'S0095',... o14
    'S0090',... o15
    };

% UCI subject IDs
subjectIDs_UCI = {'045',... o1
    '047',... o2
    '048',... o3
    '049',... o4
    '050',... o5
    '051',... o6
    '052',... o7
    '054',... o8
    '056',... o9
    '057',... o10
    '058',... o11
    '059',... o12
    '060',... o13
    '061',... o14
    '062',... o15
    }; 

subjectIDs = [subjectIDs_BU subjectIDs_UCI];

%% Setup
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'png'; % 'svg 

s.exptName      = 'cue gabor detection'; % 'cued texture detection'
s.exptShortName = 'twcf_cue_gab_det'; % 'twcf_cue_tex_det' 
s.exptFolder    = 'main_expt'; % 'main_expt'

s.site          = 'UCI+BU'; % 'BU' 'UCI' 'UCI+BU'
sites           = [ repmat({'BU'},[1,numel(subjectIDs_BU)]) repmat({'UCI'},[1,numel(subjectIDs_UCI)]) ]; 
s.sites         = sites; 

unlocked = 1; % different directories depending on whether locked or unlocked data repositories are used 

%% Loads and compiles session data to subject level structure 'dataAll'
% s.baseDir       = '/Users/karen/Dropbox/github/TWCF_FOHO'; % lab computer 
s.baseDir       = '/Users/kantian/Dropbox/github/TWCF_FOHO'; % laptop
if unlocked
    s.dataRootName  = 'twcf_expt1_data_unlocked'; 
else
    s.dataRootName  = 'twcf_expt1_data_BU'; 
end
% s.dataRootDir   = sprintf('%s/%s',s.baseDir); 
addpath(genpath(s.baseDir))

s.analDir       = sprintf('%s/twcf_expt1_analysis_BU/%s',s.baseDir,s.exptShortName); 
cd(s.analDir)

s.groupFigDir = sprintf('%s/figs',s.analDir);
if ~exist(s.groupFigDir, 'dir')
    mkdir(s.groupFigDir)
end

s.groupCSVDir = sprintf('%s/twcf_expt1_stats_BU/%s/data',s.baseDir,s.exptShortName);
if ~exist(s.groupCSVDir, 'dir')
    mkdir(s.groupCSVDir)
end

% Get figure and psychometric function fitting parameters 
p = twcf_analysisParams; 

unlocked = 1; 

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
clear a PSI
for i = 1:numel(subjectIDs)
   
    % if data unlocked, save to BU and UCI shared repo 
    if unlocked
        s.dataFolder    = sprintf('%s/%s/%s', s.baseDir, s.dataRootName, sites{i});
    else % if locked data only, save to private BU repo 
        s.dataFolder    = sprintf('%s/%s', s.baseDir, s.dataRootName);
    end

    subjectID = subjectIDs{i};
    plotPSI = 1; % no plotting of PSI until data unlocked 
    if plotPSI 
        [a(i), PSI(i)] = twcf_behavAnalysis(dataAll, subjectIDs, subjectID, s, plotPSI); 
    else
        [a(i)] = twcf_behavAnalysis(dataAll, subjectIDs, subjectID, s, plotPSI);
    end
    
    close all
end

%% Do group analysis! 
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'svg'; % svg
s.exportcsv     = 0; 

cFieldNames = {'C1','C2','CX'}; % saw grating, saw orientation, saw grating not orientation
clear g
for iC = 1:numel(cFieldNames)
    g(iC) = twcf_plotGroup(PSI,subjectIDs,s,cFieldNames,iC,sites,a);
end

%% save workspace here %% (if saving anything new, otherwise load .mat
% save('1.3data_subjects.mat')

% load('1.3data_subjects.mat')

%% Fit to group avg
addpath(genpath(s.baseDir))
cd(s.analDir)
s.saveFigs      = 1;
s.figType       = 'pdf'; 
for iC = 1 % :2 
    [g2] = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g);
end

%% Prepare export (single trial) for R


