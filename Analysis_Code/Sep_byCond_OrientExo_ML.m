

%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% /////////////////////////////////////////////////////////////////////////
%% Load data if saved previously
load([pwd '/data_out_all_' exp.settings '.mat']);
% 
% % load specific EEG dataset to make EEGLab happy
% EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);

% /////////////////////////////////////////////////////////////////////////
%% Load previously processed data segments
% /////////////////////////////////////////////////////////////////////////
% note that some settings for loading data is specified in the function itself
[ALLEEG,EEG] = LoadProcData_OrientExo_ML(exp);

eeglab redraw %reopen eeglab gui

% /////////////////////////////////////////////////////////////////////////


% /////////////////////////////////////////////////////////////////////////
%% Extract data from ALLEEG format
% /////////////////////////////////////////////////////////////////////////

% Baseline removal?
rm_baseline = 1; % 0 = no, 1 = yes
rm_timerange = [-1600 -1400]; % [min_ms max_ms] baseline latency range in ms. 
                              % [] input -> Use whole epoch as baseline

data_out_L = cell(1,length(exp.participants));    %pre-allocate
data_out_R = cell(1,length(exp.participants));    %pre-allocate
i_part = 0; %counter of subjects
for ii = 1:2:length(exp.participants)*2 %two sets per subj
    
    i_part = i_part + 1;%counter of subjects
    
    if rm_baseline == 1
        ALLEEG(ii) = pop_rmbase(ALLEEG(ii), rm_timerange);
        ALLEEG(ii+1) = pop_rmbase(ALLEEG(ii+1), rm_timerange);
    end
    
    ALLEEG(1) = pop_rmbase(ALLEEG(1),[-1600 -1400]);
    
    % Extract segments by target side
    data_out_L{i_part}(:,:,:) = squeeze(ALLEEG(ii).data(:,:,:));
    data_out_R{i_part}(:,:,:) = squeeze(ALLEEG(ii+1).data(:,:,:)); %2nd set
        
end
clear i_part ii rm_timerange


% Set-up and save relevant information
times_out = ALLEEG(1).times; % Get times variable to save
chanlocs = ALLEEG(1).chanlocs; % Get times variable to save
s_rate = ALLEEG(1).srate; % Get times variable to save

if rm_baseline == 1
    namesave = '/data_out_side_rmBL_';
else
    namesave = '/data_out_side_';
end

%save with version for large files
save([pwd namesave exp.settings '.mat'],'data_out_L','data_out_R',...
    'times_out','chanlocs','s_rate','-v7.3')
clear namesave


% /////////////////////////////////////////////////////////////////////////
%% Separate by target and cue
% /////////////////////////////////////////////////////////////////////////

data_out_Liv = cell(1,length(exp.participants));    %pre-allocate
data_out_Riv = cell(1,length(exp.participants));    %pre-allocate
data_out_Lv = cell(1,length(exp.participants));    %pre-allocate
data_out_Rv = cell(1,length(exp.participants));    %pre-allocate
for i_part = 1:length(exp.participants)
    
    % extract participant data
    tmp_L = data_out_L{i_part};
    tmp_R = data_out_R{i_part};
    
    % By left target & cue type
    data_out_Lv{i_part} = squeeze(tmp_L(:,:,(valid{i_part}(position{i_part}==0)==1)));
    data_out_Liv{i_part} = squeeze(tmp_L(:,:,(valid{i_part}(position{i_part}==0)==0)));
    % By right target & cue type
    data_out_Rv{i_part} = squeeze(tmp_R(:,:,(valid{i_part}(position{i_part}==1)==1)));
    data_out_Riv{i_part} = squeeze(tmp_R(:,:,(valid{i_part}(position{i_part}==1)==0)));
    
    clear tmp_L tmp_R    
end
clear i_part


% Set-up and save relevant information
times_out = ALLEEG(1).times; % Get times variable to save
chanlocs = ALLEEG(1).chanlocs; % Get times variable to save
s_rate = ALLEEG(1).srate; % Get times variable to save

if rm_baseline == 1
    namesave = '/data_out_cond_rmBL_';
else
    namesave = '/data_out_cond_';
end
%save with version for large files
save([pwd namesave exp.settings '.mat'],'data_out_Lv','data_out_Liv',...
    'data_out_Rv','data_out_Riv','times_out','chanlocs','s_rate','-v7.3')
clear namesave rm_baseline


