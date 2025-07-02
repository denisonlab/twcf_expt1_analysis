function [dataAll,subjectIDs] = behav_compile(subjectIDs, s)
% twcf fmri experiments
% Compiles behavioral data sessions to subjects and group data structure 

%% Load all data file names 

% * * * change here * * *
% dataFolder = "../../twcf_expt1_data_BU/twcf_cue_tex_dis/";
addpath(genpath(pwd))
% * * * * * *

%%
exptName      = s.exptName; % 'cue texture detection';
exptShortName = s.exptShortName; % 'twcf_cue_tex_det'; 
exptFolder    = s.exptFolder; % 'main_expt';

dataAll = []; 
for i = 1:numel(subjectIDs)
    switch s.exptShortName
        case {'twcf_cue_tex_det','twcf_cue_tex_dis','twcf_cue_gab_dis',...
                'twcf_cue_tex_det_fmri_pilot2'}
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
        case 'twcf_cue_tex_det_fmri'
            s.dataFolder = sprintf('%s/twcf_expt1_data_BU/%s',s.baseDir);
            subjectFolder{i} = sprintf('%s/%s/%s/%s/',s.dataFolder,s.exptShortName,subjectIDs{i},s.exptFolder); 
    end
    cd(s.dataFolder)
    addpath(genpath(pwd))

    cd(subjectFolder{i})
    listing = dir;
    dataAll(i).subjectID = subjectIDs{i}; 
    count = 1; 
    for iFile = 1:numel(listing)
        if strcmp(s.exptShortName,'twcf_cue_tex_det_fmri') && strcmp(s.exptFolder,'main_expt') % for fmri, look for data with 5 (1 session) or 10 (2 session) blocks
            if contains(listing(iFile).name,'.mat') && contains(listing(iFile).name,subjectIDs{i})
            % incomplete data
            % if contains(listing(iFile).name,'.mat') && (contains(listing(iFile).name,'trial_81_of_block_10') || contains(listing(iFile).name,'trial_81_of_block_5')) && contains(listing(iFile).name,subjectIDs{i})
                disp(['valid file ' listing(iFile).name])
                dataAll(i).dataFileNames{count} = listing(iFile).name;
                count = count + 1;
            else
                disp(['invalid file ' listing(iFile).name])
            end
        elseif strcmp(s.exptShortName,'twcf_cue_tex_det_fmri') && strcmp(s.exptFolder,'validation')
            if contains(listing(iFile).name,'.mat') && ~contains(listing(iFile).name,'trial') && contains(listing(iFile).name,subjectIDs{i})
                disp(['valid file ' listing(iFile).name])
                dataAll(i).dataFileNames{count} = listing(iFile).name;
                count = count + 1;
            else
                disp(['invalid file ' listing(iFile).name])
            end
        else % for behavior, look for completed data
            if contains(listing(iFile).name,'.mat') && ~contains(listing(iFile).name,'trial') && contains(listing(iFile).name,subjectIDs{i})
                disp(['valid file ' listing(iFile).name])
                dataAll(i).dataFileNames{count} = listing(iFile).name;
                count = count + 1;
            else
                disp(['invalid file ' listing(iFile).name])
            end
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
        
        subjectID  = dataFileSplitName{6};  % eg 'S0006';
        date       = dataFileSplitName{9};  % eg '20220614_153756';
        time       = dataFileSplitName{10}; % eg '20220614_153756';
        block      = sscanf(dataFileSplitName{end},'%d.mat'); 

        load(sprintf('%s/%s',subjectFolder{i},dataFile))
      
        nTrials = numel(data.i_trial); 
        data.sessions = ones(1,nTrials)*iFile; 
        dataAll(i).BGcolor(iFile) = p.stim.BGcolor; 
        dataAll(i).nBlocks(iFile) = block; 
        
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

    end
end

