function plotThreshTiming(dataAllThresh,s)
% Compiles duration of threshold block across observers 

% Check inputs 
if nargin<2
    s.plotFigs = 0; 
    s.saveFigs = 0; 
    s.exptName = []; 
end

%% Make figure directory 
if s.saveFigs 
    figFolder = sprintf('%s/analysis',s.groupFigDir);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder)
    end
end

%% Calculate time elapsed from first to last trial 
for iS = 1:numel(dataAllThresh)
    t1 = dataAllThresh(iS).timing.postcue(end).VBL(1);
    tend = dataAllThresh(iS).timing.postcue(end).VBL(end);
    tDur(iS) = tend-t1; % seconds 
end

%% Plot histogram of duration (in min) 
if s.plotFigs
    figure
    figureStyle
    histogram(tDur/60)
    ylabel('Count')
    xlabel('Duration of threshold (min)')
    title(s.exptName)

    if s.saveFigs
        figTitle = sprintf('%s_threshold_duration_histogram',s.exptShortName);
        saveas(gcf,sprintf('%s/%s.%s', figFolder, figTitle, s.figType))
    end
end
