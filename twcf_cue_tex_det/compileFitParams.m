% Compile psychometric function fit params
% summary statistics for Brian

fieldNames = {'detectionSawFigure','detectionSawShape','discriminationpCorrect'};
for iF = 1:numel(fieldNames)
    for iS = 1:numel(a)
        fits.(fieldNames{iF}).params(iS,:,:,:) = a(iS).(fieldNames{iF}).paramsValues;
    end
    params = fits.(fieldNames{iF}).params; 
    params(isinf(params)) = NaN; 
    fit.(fieldNames{iF}).fittedParamsMean = mean(params,1,'omitnan'); % average subjects
    fit.(fieldNames{iF}).fittedParamsMax = max(params,[],1,'omitnan'); % max 
    fit.(fieldNames{iF}).fittedParamsMin = min(params,[],1,'omitnan');
    fit.(fieldNames{iF}).fittedParamsMedian = median(params,1,'omitnan');
    
    % threshold, slope, guess-rate, lapse-rate
    % Settings
    fit.(fieldNames{iF}).PF = a(1).(fieldNames{iF}).PF; 
    fit.(fieldNames{iF}).searchGrid = a(1).(fieldNames{iF}).searchGrid; 
    fit.(fieldNames{iF}).paramValueNames = a(1).(fieldNames{iF}).paramValueNames; 
    fit.(fieldNames{iF}).paramsFree = a(1).(fieldNames{iF}).paramsFree; 
end

%% Save summary mat 
save('twcf_cue_tex_det_fitSummary.mat','fit')
