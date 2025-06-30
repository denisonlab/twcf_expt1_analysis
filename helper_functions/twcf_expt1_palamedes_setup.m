% code for setting up psychometric function fitting for twcf_expt1

switch taskType
    case 'detection',      searchGrid.gamma = 0;
    case 'discrimination', searchGrid.gamma = 0.5;
end

searchGrid.alpha    = .05:.05:3;
searchGrid.beta     = 10.^[-1:.1:1];
searchGrid.lambda   = 0:.001:.1;

paramsFree          = [1 1 0 1];

stimLevels          = p.stim.periph.lineLength_inDeg_list;
logStimLevels       = log10(stimLevels);

searchGridlog       = searchGrid;
searchGridlog.alpha = log10(searchGridlog.alpha);

[paramsFitted, logL] = PAL_PFML_Fit(logStimLevels, numPos, outOfNum, searchGridlog, paramsFree, @PAL_Gumbel);