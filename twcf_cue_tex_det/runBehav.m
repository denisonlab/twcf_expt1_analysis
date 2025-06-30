function runBehav
% Behavioral analysis for TWCF Expt 1.1, cue_tex_det 

% Finds and combines session data to subject level 
% And subjects into group structure 

% BU subject IDs 
subjectIDs_BU = {'S0002',...
    'S0004',...
    'S0005',...
    'S0006',...
    'S0007',...
    'S0008',...
    'S0010',...
    'S0011',...
    'S0012',...
    'S0013',...
    'S0014',...
    'S0015',...
    'S0016',...
    'S0017',...
    'S0018'};

% UCI subject IDs
subjectIDs_UCI = {'001',...
    '002',...
    '005',...
    '006',...
    '017',...
    '019',...
    '020',...
    '021',...
    '022',...
    '023',...
    '024',...
    '025',...
    '026',...
    '027',...
    '039'}; 

subjectIDs = [subjectIDs_BU subjectIDs_UCI];

%% Setup
s.plotFigs      = 0; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 0; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'png'; % 'svg 

s.exptName      = 'cue texture detection'; % 'cued texture detection'
s.exptShortName = 'twcf_cue_tex_det'; % 'twcf_cue_tex_det' 
s.exptFolder    = 'main_expt'; % 'main_expt'

s.site          = 'UCI+BU'; % 'BU' 'UCI'
sites           = [ repmat({'BU'},[1,numel(subjectIDs_BU)]) repmat({'UCI'},[1,numel(subjectIDs_UCI)]) ];  
s.sites         = sites; 

unlocked = 1; % different directories depending on whether locked or unlocked data repositories are used 

%% Loads and compiles session data to subject level structure 'dataAll'
% s.baseDir       = '/Users/karen/Dropbox/github/TWCF_FOHO'; % lab computer 
s.baseDir       = '/Users/kantian/Dropbox/github/TWCF_FOHO'; % laptop
if unlocked % unlocked UCI and BU data directory 
    s.dataRootName  = 'twcf_expt1_data_unlocked'; 
else % locked BU directory only 
    s.dataRootName  = 'twcf_expt1_data_BU'; 
end

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

%% Compiles data 
[dataAll,subjectIDs] = behav_compile(subjectIDs, s); % sites 
cd(s.analDir)

%% Find bg 127 (old) vs 217 (new)
% 7 subjects run with old, 23 with new
for iS=1:30
    bg(iS) = dataAll(iS).BGcolor(1); 
end
[cnt_bg, unique_bg] = hist(bg,unique(bg));

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
    s.dataFolder    = sprintf('%s/twcf_expt1_data_unlocked/%s',s.baseDir,sites{i});
    subjectID = subjectIDs{i};
    [a(i),PSI(i)] = twcf_behavAnalysis(dataAll,subjectIDs,subjectID,s); % PSI(i)
    close all % close figs 
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
% save('1.1data_subjects.mat')

% load('1.1data_subjects.mat')

%% Group figs
addpath(genpath(s.baseDir))
cd(s.analDir)
s.plotFigs      = 1; % if 1, will plot figures, otherwise just returns analysis structure 
s.saveFigs      = 1; % if 1,  will save figures to subject data directory, otherwise plots without saving  
s.figType       = 'pdf'; % svg
for iC = 1 % 1:2 % :3 % plot just C as "saw fig" 
    g2 = twcf_plotGroupAvg(a,PSI,subjectIDs,s,cFieldNames,iC,sites,g);
end



