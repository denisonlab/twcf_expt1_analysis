% Compile psychometric function fit params
% summary statistics for Brian

% Load 1.3data_subjects.mat

%% 
fieldNames = {'C1','C2'}; % 'discriminationpCorrect'
% C1 is saw figure; C2 is saw shape 
for iF = 1:numel(fieldNames)
    for iS = 1:numel(a)
        fits.(fieldNames{iF}).params(iS,:,:) = a(iS).(fieldNames{iF}).paramsFitted;
    end
    
    params = fits.(fieldNames{iF}).params; 
    params(isinf(params)) = NaN; 
    
    fits.(fieldNames{iF}).fittedParamsMean = mean(params,1,'omitnan'); % average subjects
    fits.(fieldNames{iF}).fittedParamsMax = max(params,[],1,'omitnan'); % max 
    fits.(fieldNames{iF}).fittedParamsMin = min(params,[],1,'omitnan');
    fits.(fieldNames{iF}).fittedParamsMedian = median(params,1,'omitnan');
    
    % threshold, slope, guess-rate, lapse-rate
    % Settings
    fits.(fieldNames{iF}).PF = a(1).(fieldNames{iF}).PF; 
    fits.(fieldNames{iF}).searchGrid = a(1).(fieldNames{iF}).searchGrid; 
    fits.(fieldNames{iF}).paramValueNames = a(1).(fieldNames{iF}).paramValueNames; 
    fits.(fieldNames{iF}).paramsFree = a(1).(fieldNames{iF}).paramsFree; 
end

%% Save summary mat 
save('twcf_cue_gab_det_fitSummary.mat','fits')
