function twcf_runBehav
% Behavioral analysis for TWCF Expt 1.2, cue_tex_dis

% Finds and combines session data to subject level 
% And subjects into group structure 

% BU subject IDs (n = 14) 
subjectIDs_BU = {'S0010',... % o1
    'S0021',... % o2
    'S0022',... % o3
    'S0023',... % o4
    'S0025',... % o5
    'S0026',... % o6
    'S0035',... % o7
    'S0049',... % o8
    'S0051',... % o9
    'S0054',... % o10
    'S0055',... % o11
    'S0056',... % o12
    'S0064',... % o14
    'S0065',... % o15
    'S0057', % o13 excluded for at chance A(x) performance, loading just for threshold data 
    };

%     'S0057',... % o13 excluded for at chance A(x) performance 
% S0035 has issues with fit.. need to constrain 

%     'S0024',... % XXX
%     'S0027',... % XXX
%     'S0059',... % XXX
%     'S0041',... % X
%     'S0042',... % X
%     'S0043',... % X
%     'S0060',... % X

% UCI subject IDs (n = 15) 
subjectIDs_UCI = {'027',... % o1
    '029',... % o2
    '030',... % o3
    '031',... % o4
    '032',... % o5
    '033',... % o6
    '034',... % o7
    '035',... % o8
    '037',... % o9
    '038',... % o10
    '040',... % o11
    '041',... % o12 
    '042',... % o13 
    '043',... % o14 
    '044'}; % o15

subjectIDs = [subjectIDs_BU subjectIDs_UCI];

%% Setup
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'png'; % 'svg 

s.exptName      = 'cue texture discrimination'; % 'cued texture detection'
s.exptShortName = 'twcf_cue_tex_dis'; % 'twcf_cue_tex_det' 
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

%% Compiles data 
clear dataAll
[dataAll,subjectIDs] = behav_compile(subjectIDs, s);
cd(s.analDir)

%% Compile threshold data
s.exptFolder = 'thresholding'; 
[dataAllThresh,subjectIDs] = behav_compileThresh(subjectIDs, s); % sites 
cd(s.analDir)

s.plotFigs = 0;
s.saveFigs = 0; 
s.figType = 'png'; 
plotThreshTiming(dataAllThresh,s)

%% Plot discrim detect by contrast, fit curve 
% % * * * change here * * *
% subjectID = 'S0004'; % subject of interest
% % * * * * * * * * * * * *
clear a; clear PSI;
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
s.figType       = 'png'; % 
s.exportcsv     = 1; 
cFieldNames     = {'C'}; 
clear g 
for iC = 1:numel(cFieldNames)
    g(iC) = twcf_plotGroup(PSI,subjectIDs,s,cFieldNames,iC,sites,a);
end

%% save workspace here %% 
% save('1.2data_subjects.mat')

% load('1.2data_subjects.mat')

%% Fit to group avg 
addpath(genpath(s.baseDir))
cd(s.analDir)
iC              = 1; 
s.saveFigs      = 1;
s.figType       = 'pdf'; 
[g2] = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g);


