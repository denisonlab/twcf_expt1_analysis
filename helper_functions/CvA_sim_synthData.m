function data = CvA_sim_synthData(params_A, params_C, stimLevels, nTrialsPerStim)

%% get psychometric function values at each stim level

pCorr = PAL_Weibull(params_A, stimLevels);
pVis  = PAL_Weibull(params_C, stimLevels);

% figure; hold on
% plot(stimLevels, pCorr, 'go-')
% plot(stimLevels, pVis, 'g*-')

%% make synthetic data

data.stimLevel = [];
data.acc       = [];
data.vis       = [];
for i_stim = 1:length(stimLevels)
    
    % simulate data at each stim level using psychometric function
    acc = rand(1, nTrialsPerStim) < pCorr(i_stim);
    vis = rand(1, nTrialsPerStim) < pVis(i_stim);
    
    data.stimLevel = [data.stimLevel, stimLevels(i_stim)*ones(1, nTrialsPerStim)];
    data.acc       = [data.acc, acc];
    data.vis       = [data.vis, vis];
    
    % compute outOfNum and numPos
    data.outOfNum(i_stim) = nTrialsPerStim;
    data.numPos_A(i_stim) = sum(acc);
    data.numPos_C(i_stim) = sum(vis);
end

