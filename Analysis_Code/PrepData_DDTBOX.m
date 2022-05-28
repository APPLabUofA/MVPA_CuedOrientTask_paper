
% Clears the workspace and closes all figure windows
ccc 

% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);
load([pwd '/BEHerr_M_' exp.settings '.mat']);

% load specific EEG dataset to make EEGLab happy
% EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Load extracted features
% load('feats_times876_byTargets_ML_v1.mat')
load('feats_v2_times751_byTargets_ML_v1.mat')

% /////////////////////////////////////////////////////////////////////////
%% Folder to save the subject data

% bdir = [pwd '\decode_v1\' exp.settings '\'];
bdir = [pwd '\decode_v2\' exp.settings '\'];
% if folder doesn't exist yet, create one
if ~exist(bdir)
    mkdir(bdir);
end

% Folder for AAC
bdir_AAC = [bdir 'AAC\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_AAC)
    mkdir(bdir_AAC);
end

% Folder for lowFreq_amp
bdir_lowFreq_amp = [bdir 'lowFreq_amp\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_lowFreq_amp)
    mkdir(bdir_lowFreq_amp);
end

% Folder for PACcan
bdir_PACcan = [bdir 'PACcan\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PACcan)
    mkdir(bdir_PACcan);
end

% Folder for PACoz
bdir_PACoz = [bdir 'PACoz\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PACoz)
    mkdir(bdir_PACoz);
end

% Folder for PACplv
bdir_PACplv = [bdir 'PACplv\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PACplv)
    mkdir(bdir_PACplv);
end

% Folder for PACtort
bdir_PACtort = [bdir 'PACtort\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PACtort)
    mkdir(bdir_PACtort);
end

% Folder for PACtort
bdir_PACglm = [bdir 'PACglm\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PACglm)
    mkdir(bdir_PACglm);
end



% /////////////////////////////////////////////////////////////////////////
%% Save channel information
chaninfo = EEG.chaninfo;
chanlocs = EEG.chanlocs;
save([bdir 'channel_OrientExo_ML'],'chaninfo','chanlocs')


% /////////////////////////////////////////////////////////////////////////
%% Extract features and save for DDTBOX

% Create info_condition cell which will be saved with each subj
info_condition = {}; %rest
info_condition(1,:) = {'Condition:',1,2,3,4};
info_condition(2,:) = {'Type:','L_NI','L_I','R_NI','R_I'};

% Data needs to be organized by participant:
% eeg_sorted_cond{run, condition}(timepoints, channels, epochs)
% SVR_labels{run, condition}(epoch_number) named sbj#_regress_sorted_data

for ipart = 1:length(exp.participants)
    
    %% Support Vector Regression (SVR) condition labels
    SVR_labels = {}; % reset data cell
    
%     SVR_labels{1,1} = squeeze(errordeg{ipart}(valid{ipart}(position{ipart}==0)==0))'; %Liv
%     SVR_labels{1,2} = squeeze(errordeg{ipart}(valid{ipart}(position{ipart}==0)==1))'; %Lv
%     SVR_labels{1,3} = squeeze(errordeg{ipart}(valid{ipart}(position{ipart}==1)==0))'; %Riv
%     SVR_labels{1,4} = squeeze(errordeg{ipart}(valid{ipart}(position{ipart}==1)==1))'; %Rv
    
    % Absolute errors
    SVR_labels{1,1} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==0)==0)))'; %Liv
    SVR_labels{1,2} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==0)==1)))'; %Lv
    SVR_labels{1,3} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==1)==0)))'; %Riv
    SVR_labels{1,4} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==1)==1)))'; %Rv
    
      % Probability trial is from the guess distribution 
%     SVR_labels{1,1} = squeeze(M_invalidL{ipart}(:,2)); %Liv
%     SVR_labels{1,2} = squeeze(M_validL{ipart}(:,2)); %Lv
%     SVR_labels{1,3} = squeeze(M_invalidR{ipart}(:,2)); %Riv
%     SVR_labels{1,4} = squeeze(M_validR{ipart}(:,2)); %Rv
    
    
    %% AAC
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_AAC_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_AAC_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_AAC_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_AAC_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_AAC 'EEG_data'])
        mkdir([bdir_AAC 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_AAC 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_AAC 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_AAC 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% lowFreq amp
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data (no reordering)
%     eeg_sorted_cond(1).data(:,:,:) = feat_lowFreq_amp_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_lowFreq_amp_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_lowFreq_amp_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_lowFreq_amp_Rv{ipart,1}; 
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_lowFreq_amp 'EEG_data'])
        mkdir([bdir_lowFreq_amp 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_lowFreq_amp 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_lowFreq_amp 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_lowFreq_amp 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PACcan
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_PACcan_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_PACcan_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_PACcan_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_PACcan_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PACcan 'EEG_data'])
        mkdir([bdir_PACcan 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PACcan 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_PACcan 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_PACcan 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PACoz
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_PACoz_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_PACoz_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_PACoz_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_PACoz_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PACoz 'EEG_data'])
        mkdir([bdir_PACoz 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PACoz 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_PACoz 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_PACoz 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PACplv
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_PACplv_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_PACplv_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_PACplv_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_PACplv_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PACplv 'EEG_data'])
        mkdir([bdir_PACplv 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PACplv 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_PACplv 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_PACplv 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PACtort
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_PACtort_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_PACtort_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_PACtort_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_PACtort_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PACtort 'EEG_data'])
        mkdir([bdir_PACtort 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PACtort 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_PACtort 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_PACtort 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
     %% PACglm
%     eeg_sorted_cond = {}; % reset data cell
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = permute(feat_PACglm_Liv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(2).data(:,:,:) = permute(feat_PACglm_Lv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(3).data(:,:,:) = permute(feat_PACglm_Riv{ipart,1},[3 2 1]); %reorder dimensions
%     eeg_sorted_cond(4).data(:,:,:) = permute(feat_PACglm_Rv{ipart,1},[3 2 1]); %reorder dimensions
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PACglm 'EEG_data'])
        mkdir([bdir_PACglm 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PACglm 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
    save([bdir_PACglm 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
%     save([bdir_PACglm 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
end
clear ipart eeg_sorted_cond info_condition SVR_labels
































