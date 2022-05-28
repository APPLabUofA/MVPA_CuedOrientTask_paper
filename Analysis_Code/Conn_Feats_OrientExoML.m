%% Clean up workspace
ccc

%% Add necessary folders to path list
addpath 'C:\Users\ssshe\Documents\MathLab\Analysis\OrientWheel_Exo' %specific to my computer
%location of fieldtrip functions
addpath 'C:\Users\ssshe\Documents\MathLab\Experiments\matlab\fieldtrip-20160928' %specific to my computer
addpath 'C:\Users\ssshe\Documents\MathLab\Experiments\matlab\fieldtrip-20160928\preproc' %specific to my computer

% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% Change exp structure name to exp2 cuz exp function
exp2 = exp;
clear exp

% /////////////////////////////////////////////////////////////////////////
%% Load EEG data
load([pwd '/data_out_cond_' exp2.settings '.mat']);

% load specific EEG dataset to make EEGLab happy
EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Choose frequency bands for analysis
freqband = [8,14]; %alpha

% /////////////////////////////////////////////////////////////////////////
%% Set-up time parameters for analysis
% Define time period in ms (relative to aligned event)
% timewindow = [ 200 1200 ]; % cue-aligned
timewindow = [ -1000 0 ]; % target-aligned
winSize = (1/freqband(1))/(1/s_rate)*3; % must be at least 1/lowFreq < window_length
winStep = 5; %size of steps time window moves across trial, ms (e.g., 50)

% /////////////////////////////////////////////////////////////////////////
%% Set other parameters

% Choose electrodes
select_elect = exp2.brainelecs; %all brain electrodes

% For the ft_preproc_bandpassfilter
filt_order = 4; %default

Fs = s_rate; %sampling rate
times = times_out; %time points of epoch

% Add some padding to time to keep dimesions the same
times = ft_preproc_padding(times, 'nan', 600);


% /////////////////////////////////////////////////////////////////////////
%% Laplacian transformation (CSD)
% /////////////////////////////////////////////////////////////////////////

% extract XYZ coordinates from EEG structure
X = [EEG.chanlocs(exp2.brainelecs).X];
Y = [EEG.chanlocs(exp2.brainelecs).Y];
Z = [EEG.chanlocs(exp2.brainelecs).Z];

% surf_lapN = laplacian_nola(X,Y,Z,EEGsig,100);
% surf_lapP = laplacian_perrinX(EEGsig,X,Y,Z,[],1e-5); 

lapP_Liv = cellfun(@(sig) laplacian_perrinX(sig(exp2.brainelecs,:,:),X,Y,Z,[],1e-5),data_out_Liv,'UniformOutput',false);
lapP_Lv = cellfun(@(sig) laplacian_perrinX(sig(exp2.brainelecs,:,:),X,Y,Z,[],1e-5),data_out_Lv,'UniformOutput',false);
lapP_Riv = cellfun(@(sig) laplacian_perrinX(sig(exp2.brainelecs,:,:),X,Y,Z,[],1e-5),data_out_Riv,'UniformOutput',false);
lapP_Rv = cellfun(@(sig) laplacian_perrinX(sig(exp2.brainelecs,:,:),X,Y,Z,[],1e-5),data_out_Rv,'UniformOutput',false);

clear X Y Z

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Connectivity measures
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////

connfeat_plv_Liv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_corr_Liv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_pli_Liv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_wpli_Liv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_dwpli_Liv  = cell(length(exp2.participants),1); %pre-allocate
connfeat_icoh_Liv   = cell(length(exp2.participants),1); %pre-allocate

connfeat_plv_Lv     = cell(length(exp2.participants),1); %pre-allocate
connfeat_corr_Lv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_pli_Lv     = cell(length(exp2.participants),1); %pre-allocate
connfeat_wpli_Lv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_dwpli_Lv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_icoh_Lv    = cell(length(exp2.participants),1); %pre-allocate

connfeat_plv_Riv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_corr_Riv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_pli_Riv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_wpli_Riv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_dwpli_Riv  = cell(length(exp2.participants),1); %pre-allocate
connfeat_icoh_Riv   = cell(length(exp2.participants),1); %pre-allocate

connfeat_plv_Rv     = cell(length(exp2.participants),1); %pre-allocate
connfeat_corr_Rv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_pli_Rv     = cell(length(exp2.participants),1); %pre-allocate
connfeat_wpli_Rv    = cell(length(exp2.participants),1); %pre-allocate
connfeat_dwpli_Rv   = cell(length(exp2.participants),1); %pre-allocate
connfeat_icoh_Rv    = cell(length(exp2.participants),1); %pre-allocate

for i_part = 1:length(exp2.participants)

    % Select single subject data   
    partdata_Liv = data_out_Liv{i_part}(exp2.brainelecs,:,:); % chan x time x trial
    partdata_Lv = data_out_Lv{i_part}(exp2.brainelecs,:,:); % chan x time x trial
    partdata_Riv = data_out_Riv{i_part}(exp2.brainelecs,:,:); % chan x time x trial
    partdata_Rv = data_out_Rv{i_part}(exp2.brainelecs,:,:); % chan x time x trial
    
    % Select single subject data - Laplacian transformed   
    partlapP_Liv = lapP_Liv{i_part}; % chan x time x trial
    partlapP_Lv = lapP_Lv{i_part}; % chan x time x trial
    partlapP_Riv = lapP_Riv{i_part}; % chan x time x trial
    partlapP_Rv = lapP_Rv{i_part}; % chan x time x trial
    
    
    %% Bandpass Filter
    filtdata_Liv = NaN(length(exp2.brainelecs),size(partdata_Liv,3),length(times)); %pre-allocate
    filtdata_Lv = NaN(length(exp2.brainelecs),size(partdata_Lv,3),length(times)); %pre-allocate
    filtdata_Riv = NaN(length(exp2.brainelecs),size(partdata_Riv,3),length(times)); %pre-allocate
    filtdata_Rv = NaN(length(exp2.brainelecs),size(partdata_Rv,3),length(times)); %pre-allocate
    filtlapP_Liv = NaN(length(exp2.brainelecs),size(partdata_Liv,3),length(times)); %pre-allocate
    filtlapP_Lv = NaN(length(exp2.brainelecs),size(partdata_Lv,3),length(times)); %pre-allocate
    filtlapP_Riv = NaN(length(exp2.brainelecs),size(partdata_Riv,3),length(times)); %pre-allocate
    filtlapP_Rv = NaN(length(exp2.brainelecs),size(partdata_Rv,3),length(times)); %pre-allocate
    for ichan = 1:length(exp2.brainelecs)
        
        % Add some padding before and after data
        singleChanEEG_Liv = ft_preproc_padding(squeeze(partdata_Liv(ichan,:,:))', 'zero', 600);
        singleChanEEG_Lv = ft_preproc_padding(squeeze(partdata_Lv(ichan,:,:))', 'zero', 600);
        singleChanEEG_Riv = ft_preproc_padding(squeeze(partdata_Riv(ichan,:,:))', 'zero', 600);
        singleChanEEG_Rv = ft_preproc_padding(squeeze(partdata_Rv(ichan,:,:))', 'zero', 600);
        
        % Add some padding before and after data - Laplacian transformed 
        singleChanlapP_Liv = ft_preproc_padding(squeeze(partlapP_Liv(ichan,:,:))', 'zero', 600);
        singleChanlapP_Lv = ft_preproc_padding(squeeze(partlapP_Lv(ichan,:,:))', 'zero', 600);
        singleChanlapP_Riv = ft_preproc_padding(squeeze(partlapP_Riv(ichan,:,:))', 'zero', 600);
        singleChanlapP_Rv = ft_preproc_padding(squeeze(partlapP_Rv(ichan,:,:))', 'zero', 600);
        
        
        % Bandpass filter
        filtdata_Liv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanEEG_Liv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtdata_Lv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanEEG_Lv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtdata_Riv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanEEG_Riv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtdata_Rv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanEEG_Rv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        clear singleChanEEG_*
        
        % Bandpass filter - Laplacian transformed
        filtlapP_Liv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanlapP_Liv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtlapP_Lv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanlapP_Lv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtlapP_Riv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanlapP_Riv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        filtlapP_Rv(ichan,:,:) = ft_preproc_bandpassfilter(singleChanlapP_Rv,Fs,...
            [freqband(1) freqband(2)],filt_order,'but','twopass','no');
        clear singleChanlapP_*
        
    end
    clear ichan
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set time window(s) of analysis
    % convert time window to logical
    timelim = times>=timewindow(1) & times<=timewindow(2);
    Ntime = length(times(timelim)); %length of analysis window
    
    time_idx = NaN(length(1:(Ntime-winSize)/winStep),1); %pre-allocate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ij = 1:floor((Ntime-winSize)/winStep)
                
        % Get current time points for analysis
        tt = [(ij-1)*winStep+1:(ij-1)*winStep+winSize];

        % Mid-window for time indexing
        time_idx(ij) = tt(round(winSize/2));
        
        
        % Liv -------------------------------------------------------------
        % measures using the Laplacian data
        eeglapP_tmp = filtlapP_Liv(:,:,tt);
        
        connfeat_plv_Liv{i_part}(ij,:,:,:) = conn_PLV(eeglapP_tmp);
        connfeat_corr_Liv{i_part}(ij,:,:,:) = conn_corr_amp(eeglapP_tmp);
        
        % measures using the original data
        eeg_tmp = filtdata_Liv(:,:,tt);
        
        connfeat_pli_Liv{i_part}(ij,:,:,:) = conn_PLI(eeg_tmp);
        [connfeat_wpli_Liv{i_part}(ij,:,:,:),connfeat_dwpli_Liv{i_part}(ij,:,:,:)]...
            = conn_wPLI(eeg_tmp);
        connfeat_icoh_Liv{i_part}(ij,:,:,:) = conn_iCOH(eeg_tmp);
        
        clear eeglapP_tmp eeg_tmp
        
        
        % Lv -------------------------------------------------------------
        % measures using the Laplacian data
        eeglapP_tmp = filtlapP_Lv(:,:,tt);
        
        connfeat_plv_Lv{i_part}(ij,:,:,:) = conn_PLV(eeglapP_tmp);
        connfeat_corr_Lv{i_part}(ij,:,:,:) = conn_corr_amp(eeglapP_tmp);
        
        % measures using the original data
        eeg_tmp = filtdata_Lv(:,:,tt);
        
        connfeat_pli_Lv{i_part}(ij,:,:,:) = conn_PLI(eeg_tmp);
        [connfeat_wpli_Lv{i_part}(ij,:,:,:),connfeat_dwpli_Lv{i_part}(ij,:,:,:)]...
            = conn_wPLI(eeg_tmp);
        connfeat_icoh_Lv{i_part}(ij,:,:,:) = conn_iCOH(eeg_tmp);
        
        clear eeglapP_tmp eeg_tmp
        
        
        % Riv -------------------------------------------------------------
        % measures using the Laplacian data
        eeglapP_tmp = filtlapP_Riv(:,:,tt);
        
        connfeat_plv_Riv{i_part}(ij,:,:,:) = conn_PLV(eeglapP_tmp);
        connfeat_corr_Riv{i_part}(ij,:,:,:) = conn_corr_amp(eeglapP_tmp);
        
        % measures using the original data
        eeg_tmp = filtdata_Riv(:,:,tt);
        
        connfeat_pli_Riv{i_part}(ij,:,:,:) = conn_PLI(eeg_tmp);
        [connfeat_wpli_Riv{i_part}(ij,:,:,:),connfeat_dwpli_Riv{i_part}(ij,:,:,:)]...
            = conn_wPLI(eeg_tmp);
        connfeat_icoh_Riv{i_part}(ij,:,:,:) = conn_iCOH(eeg_tmp);
        
        clear eeglapP_tmp eeg_tmp
        
        
        % Rv -------------------------------------------------------------
        % measures using the Laplacian data
        eeglapP_tmp = filtlapP_Rv(:,:,tt);
        
        connfeat_plv_Rv{i_part}(ij,:,:,:) = conn_PLV(eeglapP_tmp);
        connfeat_corr_Rv{i_part}(ij,:,:,:) = conn_corr_amp(eeglapP_tmp);
        
        % measures using the original data
        eeg_tmp = filtdata_Rv(:,:,tt);
        
        connfeat_pli_Rv{i_part}(ij,:,:,:) = conn_PLI(eeg_tmp);
        [connfeat_wpli_Rv{i_part}(ij,:,:,:),connfeat_dwpli_Rv{i_part}(ij,:,:,:)]...
            = conn_wPLI(eeg_tmp);
        connfeat_icoh_Rv{i_part}(ij,:,:,:) = conn_iCOH(eeg_tmp);
        
        clear eeglapP_tmp eeg_tmp tt 
    end
    clear ij partdata_* partlapP_* filtdata_* filtlapP_*
end
clear i_part



% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
%% Save with version for large files

tmp = times(timelim);
conn_time = tmp(time_idx); %mid-point of each data time window (for plotting)
times = times_out; %before padding

save([pwd '/connfeats_times' num2str(length(time_idx)) '_' exp2.settings '.mat'],...
    'timewindow','times','select_elect','filt_order','winSize','winStep','time_idx',...
    'conn_time','-regexp','connfeat_*','-v7.3');


% /////////////////////////////////////////////////////////////////////////
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////

clear Fs filt_order highFreq winSize winStep method Ntime data_out_* lapP_*

















