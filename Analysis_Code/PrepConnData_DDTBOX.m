
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
% load('connfeats_times125_byTargets_ML_v1.mat')

% /////////////////////////////////////////////////////////////////////////
%% Folder to save the subject data

bdir = [pwd '\conn_decode_v1\' exp.settings '\'];
% if folder doesn't exist yet, create one
if ~exist(bdir)
    mkdir(bdir);
end

% Folder for correlation
bdir_corr = [bdir 'corr\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_corr)
    mkdir(bdir_corr);
end

% Folder for dwPLI
bdir_dwPLI = [bdir 'dwPLI\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_dwPLI)
    mkdir(bdir_dwPLI);
end

% Folder for iCOH
bdir_iCOH = [bdir 'iCOH\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_iCOH)
    mkdir(bdir_iCOH);
end

% Folder for PLI
bdir_PLI = [bdir 'PLI\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PLI)
    mkdir(bdir_PLI);
end

% Folder for PLV
bdir_PLV = [bdir 'PLV\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_PLV)
    mkdir(bdir_PLV);
end

% Folder for wPLI
bdir_wPLI = [bdir 'wPLI\'];
% if folder doesn't exist yet, create one
if ~exist(bdir_wPLI)
    mkdir(bdir_wPLI);
end


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% Reorganize data

% List of electrode combos
conn_chanlocs = struct();
count = 0;
for ii = 1:(length(exp.elec_names)-1)
    for kk = (ii+1):length(exp.elec_names)
        count = 1 + count;
        conn_chanlocs(count).label{1,1} = [exp.elec_names{ii} '-' exp.elec_names{kk}];
        conn_chanlocs(count).index(1,1) = ii;
        conn_chanlocs(count).index(1,2) = kk;
    end
    clear kk
end
clear ii count


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Reorganize so electrode combos are vectors
% original: connfeat_data{participant}(time x trials x channel x channel)
% new: feat_data{participant}(time x channel combo x trials)

feat_corr_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_corr_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_corr_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_corr_Rv    = cell(length(exp.participants),1); %pre-allocate

feat_dwpli_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_dwpli_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_dwpli_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_dwpli_Rv    = cell(length(exp.participants),1); %pre-allocate

feat_icoh_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_icoh_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_icoh_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_icoh_Rv    = cell(length(exp.participants),1); %pre-allocate

feat_pli_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_pli_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_pli_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_pli_Rv    = cell(length(exp.participants),1); %pre-allocate

feat_plv_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_plv_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_plv_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_plv_Rv    = cell(length(exp.participants),1); %pre-allocate

feat_wpli_Liv   = cell(length(exp.participants),1); %pre-allocate
feat_wpli_Lv    = cell(length(exp.participants),1); %pre-allocate
feat_wpli_Riv   = cell(length(exp.participants),1); %pre-allocate
feat_wpli_Rv    = cell(length(exp.participants),1); %pre-allocate

for ipart = 1:length(exp.participants)
    count = 0;
    for ii = 1:(length(exp.elec_names)-1)
        for kk = (ii+1):length(exp.elec_names)
            count = 1 + count;
            
            feat_corr_Liv{ipart}(:,count,:) = connfeat_corr_Liv{ipart}(:,:,ii,kk);
            feat_corr_Lv{ipart}(:,count,:)  = connfeat_corr_Lv{ipart}(:,:,ii,kk);
            feat_corr_Riv{ipart}(:,count,:) = connfeat_corr_Riv{ipart}(:,:,ii,kk);
            feat_corr_Rv{ipart}(:,count,:)  = connfeat_corr_Rv{ipart}(:,:,ii,kk);
            
            feat_dwpli_Liv{ipart}(:,count,:) = connfeat_dwpli_Liv{ipart}(:,:,ii,kk);
            feat_dwpli_Lv{ipart}(:,count,:)  = connfeat_dwpli_Lv{ipart}(:,:,ii,kk);
            feat_dwpli_Riv{ipart}(:,count,:) = connfeat_dwpli_Riv{ipart}(:,:,ii,kk);
            feat_dwpli_Rv{ipart}(:,count,:)  = connfeat_dwpli_Rv{ipart}(:,:,ii,kk);
            
            feat_icoh_Liv{ipart}(:,count,:) = connfeat_icoh_Liv{ipart}(:,:,ii,kk);
            feat_icoh_Lv{ipart}(:,count,:)  = connfeat_icoh_Lv{ipart}(:,:,ii,kk);
            feat_icoh_Riv{ipart}(:,count,:) = connfeat_icoh_Riv{ipart}(:,:,ii,kk);
            feat_icoh_Rv{ipart}(:,count,:)  = connfeat_icoh_Rv{ipart}(:,:,ii,kk);
            
            feat_pli_Liv{ipart}(:,count,:) = connfeat_pli_Liv{ipart}(:,:,ii,kk);
            feat_pli_Lv{ipart}(:,count,:)  = connfeat_pli_Lv{ipart}(:,:,ii,kk);
            feat_pli_Riv{ipart}(:,count,:) = connfeat_pli_Riv{ipart}(:,:,ii,kk);
            feat_pli_Rv{ipart}(:,count,:)  = connfeat_pli_Rv{ipart}(:,:,ii,kk);
            
            feat_plv_Liv{ipart}(:,count,:) = connfeat_plv_Liv{ipart}(:,:,ii,kk);
            feat_plv_Lv{ipart}(:,count,:)  = connfeat_plv_Lv{ipart}(:,:,ii,kk);
            feat_plv_Riv{ipart}(:,count,:) = connfeat_plv_Riv{ipart}(:,:,ii,kk);
            feat_plv_Rv{ipart}(:,count,:)  = connfeat_plv_Rv{ipart}(:,:,ii,kk);
            
            feat_wpli_Liv{ipart}(:,count,:) = connfeat_wpli_Liv{ipart}(:,:,ii,kk);
            feat_wpli_Lv{ipart}(:,count,:)  = connfeat_wpli_Lv{ipart}(:,:,ii,kk);
            feat_wpli_Riv{ipart}(:,count,:) = connfeat_wpli_Riv{ipart}(:,:,ii,kk);
            feat_wpli_Rv{ipart}(:,count,:)  = connfeat_wpli_Rv{ipart}(:,:,ii,kk);
            
        end
        clear kk
    end
    clear ii count
    
    
end
clear ipart connfeat_corr_* connfeat_dwpli_* connfeat_icoh_* connfeat_pli_*...
    connfeat_plv_* connfeat_wpli_*


save([pwd '/connfeatraw_times' num2str(length(time_idx)) '_' exp2.settings '.mat'],...
    'timewindow','times','select_elect','conn_chanlocs','chanlocs','time_idx',...
    'conn_time','-regexp','feat_*','-v7.3');



% /////////////////////////////////////////////////////////////////////////
%% Save channel information
chaninfo = EEG.chaninfo;
chanlocs = EEG.chanlocs;
save([bdir 'channel_OrientExo_ML'],'chaninfo','chanlocs')
save([bdir 'connchannel_OrientExo_ML'],'conn_chanlocs')


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
%     SVR_labels{1,1} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==0)==0)))'; %Liv
%     SVR_labels{1,2} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==0)==1)))'; %Lv
%     SVR_labels{1,3} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==1)==0)))'; %Riv
%     SVR_labels{1,4} = squeeze(abs(errordeg{ipart}(valid{ipart}(position{ipart}==1)==1)))'; %Rv
    
    % Probability trial is from the guess distribution 
    SVR_labels{1,1} = squeeze(M_invalidL{ipart}(:,2)); %Liv
    SVR_labels{1,2} = squeeze(M_validL{ipart}(:,2)); %Lv
    SVR_labels{1,3} = squeeze(M_invalidR{ipart}(:,2)); %Riv
    SVR_labels{1,4} = squeeze(M_validR{ipart}(:,2)); %Rv
    
    
    %% Correlation
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_corr_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_corr_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_corr_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_corr_Rv{ipart,1}; 
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_corr 'EEG_data'])
        mkdir([bdir_corr 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_corr 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_corr 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_corr 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% dwPLI
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_dwpli_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_dwpli_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_dwpli_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_dwpli_Rv{ipart,1}; 
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_dwPLI 'EEG_data'])
        mkdir([bdir_dwPLI 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_dwPLI 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_dwPLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_dwPLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% iCOH
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_icoh_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_icoh_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_icoh_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_icoh_Rv{ipart,1};
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_iCOH 'EEG_data'])
        mkdir([bdir_iCOH 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_iCOH 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_iCOH 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_iCOH 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PLI
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_pli_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_pli_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_pli_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_pli_Rv{ipart,1};
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PLI 'EEG_data'])
        mkdir([bdir_PLI 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PLI 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_PLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_PLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% PLV
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_plv_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_plv_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_plv_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_plv_Rv{ipart,1};
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_PLV 'EEG_data'])
        mkdir([bdir_PLV 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_PLV 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_PLV 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_PLV 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
    %% wPLI
    
    % extract data
%     eeg_sorted_cond(1).data(:,:,:) = feat_wpli_Liv{ipart,1}; 
%     eeg_sorted_cond(2).data(:,:,:) = feat_wpli_Lv{ipart,1}; 
%     eeg_sorted_cond(3).data(:,:,:) = feat_wpli_Riv{ipart,1}; 
%     eeg_sorted_cond(4).data(:,:,:) = feat_wpli_Rv{ipart,1};
    
    % if folder doesn't exist yet, create one
    if ~exist([bdir_wPLI 'EEG_data'])
        mkdir([bdir_wPLI 'EEG_data']);
    end
    
    % save subj data
%     save([bdir_wPLI 'EEG_data\sbj' exp.participants{ipart}],'eeg_sorted_cond','info_condition')
%     save([bdir_wPLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_data'],'SVR_labels')
    save([bdir_wPLI 'EEG_data\sbj' exp.participants{ipart} '_regress_sorted_dataM'],'SVR_labels')
    clear eeg_sorted_cond
    
end
clear ipart eeg_sorted_cond info_condition SVR_labels
































