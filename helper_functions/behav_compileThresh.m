function [dataAll,subjectIDs] = behav_compileThresh(subjectIDs, s)
% twcf expts 1.1-1.4
% Compiles behavioral thresholding data per subject 

%% Load all data file names 

% * * * change here * * *
% dataFolder = "../../twcf_expt1_data_BU/twcf_cue_tex_dis/";
addpath(genpath(pwd))
% * * * * * *

%%
exptName      = s.exptName; % 'cued texture detection';
exptShortName = s.exptShortName; % 'twcf_cue_tex_det'; 
exptFolder    = 'thresholding'; % s.exptFolder; % 'main_expt';

dataAll = []; 
for i = 1:numel(subjectIDs)
    switch s.exptShortName
        case {'twcf_cue_tex_det','twcf_cue_tex_dis','twcf_cue_gab_dis'}
            switch s.site
                case 'BU'
                    s.dataFolder = sprintf('%s/twcf_expt1_data_BU/%s',s.baseDir);
                case 'UCI+BU'
                    s.dataFolder    = sprintf('%s/twcf_expt1_data_unlocked/%s',s.baseDir,s.sites{i});
            end
            subjectFolder{i} = sprintf('%s/%s/%s/%s/',s.dataFolder,s.exptShortName,subjectIDs{i},s.exptFolder); 
        case 'twcf_cue_gab_det'
            switch s.site 
                case 'BU'
                    s.dataFolder = sprintf('%s/twcf_expt1_data_BU/%s',s.baseDir,s.sites{i}); % locked
                    subjectFolder{i} = sprintf('%s/%s/corrected/%s/%s/',s.dataFolder,s.exptShortName,subjectIDs{i},s.exptFolder); 
                case 'UCI+BU'
                    s.dataFolder = sprintf('%s/twcf_expt1_data_unlocked/%s',s.baseDir,s.sites{i});
                    subjectFolder{i} = sprintf('%s/%s/corrected/%s/%s/',s.dataFolder,s.exptShortName,subjectIDs{i},s.exptFolder); 
            end
        case 'twcf_cue_gab_dis'
            s.dataFolder = sprintf('%s/twcf_expt1_data_BU/%s',s.baseDir);
    end
    cd(s.dataFolder)
    addpath(genpath(pwd))

    cd(subjectFolder{i})
    listing = dir;
    dataAll(i).subjectID = subjectIDs{i}; 
    count = 1; 
    for iFile = 1:numel(listing)
        if contains(listing(iFile).name,'.mat') && ~contains(listing(iFile).name,'trial') && contains(listing(iFile).name,subjectIDs{i})
            disp(['valid file ' listing(iFile).name])
            dataAll(i).dataFileNames{count} = listing(iFile).name; 
            count = count + 1; 
        else 
            disp(['invalid file ' listing(iFile).name])
        end
    end
end

%% Compile sessions per subject
for i = 1:numel(subjectIDs)
    subjectID = subjectIDs{i};
    disp(subjectID)
    for iFile = 1:numel(dataAll(i).dataFileNames)
        fprintf('file #%d',iFile)
        dataFile = dataAll(i).dataFileNames{iFile};

        dataFileSplitName = strsplit(dataFile,'_');

        % subjectID  = dataFileSplitName{6};  % eg 'S0006';
        % date       = dataFileSplitName{9};  % eg '20220614_153756';
        % time       = dataFileSplitName{10}; % eg '20220614_153756';

        load(sprintf('%s/%s',subjectFolder{i},dataFile))

        nTrials = numel(data.i_trial); 
        data.sessions = ones(1,nTrials)*iFile; 
        
        % load data 
        fields = fieldnames(data);
        for iField = 1:numel(fields) 
            if ~strcmp(fields{iField},'passCurrent') % not all sessions have passCurrent field 
                fieldData = data.(fields{iField});
                if iFile==1
                    dataAll(i).(fields{iField}) = fieldData;
                elseif iFile>1
                    catFieldData = cat(ndims(fieldData),fieldData,dataAll(i).(fields{iField}));
                    dataAll(i).(fields{iField}) = catFieldData;
                end
            end
        end

        % load timing
        % fields = fieldnames(timing);
        % for iField = 1:numel(fields) 
        %     if ~strcmp(fields{iField},'passCurrent') % not all sessions have passCurrent field 
        %         fieldData = timing.(fields{iField});
        %         if iFile==1
        %             dataAll(i).timing.(fields{iField}) = fieldData;
        %         elseif iFile>1
        %             catFieldData = cat(ndims(fieldData),fieldData,dataAll(i).timing.(fields{iField}));
        %             dataAll(i).timing.(fields{iField}) = catFieldData;
        %         end
        %     end
        % end
    end
end


