%% Clean up workspace
ccc

%% Add necessary folders to path list
addpath 'C:\Users\ssshe\Documents\MathLab\Analysis\OrientWheel_Exo' %specific to my computer
%location of different PAC functions
addpath 'C:\Users\ssshe\Documents\MathLab\Personal_Folders\Sarah\Matlab_Functions\PACmeg' %specific to my computer
%location of fieldtrip functions
addpath 'C:\Users\ssshe\Documents\MathLab\Experiments\matlab\fieldtrip-20160928' %specific to my computer
addpath 'C:\Users\ssshe\Documents\MathLab\Experiments\matlab\fieldtrip-20160928\preproc' %specific to my computer

% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

%% Load EEG data
load([pwd '/data_out_cond_' exp.settings '.mat']);

% load specific EEG dataset to make EEGLab happy
EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
eeglab redraw

%% Load oscillation info
load([pwd '/fooof_band_trialpk_' exp.settings '.mat']); %trials
load([pwd '/fooof_band_pk_' exp.settings '.mat']); %averages

%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);

% /////////////////////////////////////////////////////////////////////////
%% Choose high frequency bands for analysis

% High frequency band in Hz (e.g. [40:2:100])
    % considered amplitude frequency in PAC
% highFreq = 30:2:45; %gamma
highFreq = [30,45]; %gamma

% /////////////////////////////////////////////////////////////////////////
%% Choose PAC calculation method
% method = 'canolty';
% ('tort','ozkurt','plv','canolty')

% /////////////////////////////////////////////////////////////////////////
%% Set-up time parameters for analysis
% Define time period in ms (relative to aligned event)
% timewindow = [ 200 1200 ]; % cue-aligned
timewindow = [ -1000 0 ]; % target-aligned
winSize = 125*2; % must be 1/lowFreq < window_length
winStep = 1; %size of steps time window moves across trial, ms (e.g., 50)

% Choose electrodes
select_elect = exp.brainelecs; %all brain electrodes

% /////////////////////////////////////////////////////////////////////////
%% Set other parameters
% For the ft_preproc_bandpassfilter
filt_order = 4; %default

Fs = s_rate; %sampling rate
times = times_out; %time points of epoch

% Add some padding to time to keep dimesions the same
times = ft_preproc_padding(times, 'nan', 600);

% /////////////////////////////////////////////////////////////////////////



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through subjects. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feat_lowFreq_amp_Liv = cell(length(exp.participants),1); %pre-allocate
feat_AAC_Liv         = cell(length(exp.participants),1); %pre-allocate
feat_PACcan_Liv      = cell(length(exp.participants),1); %pre-allocate
feat_PACglm_Liv      = cell(length(exp.participants),1); %pre-allocate
feat_PACoz_Liv       = cell(length(exp.participants),1); %pre-allocate
feat_PACplv_Liv      = cell(length(exp.participants),1); %pre-allocate
feat_PACtort_Liv     = cell(length(exp.participants),1); %pre-allocate
feat_times_Liv       = cell(length(exp.participants),1); %pre-allocate

feat_lowFreq_amp_Lv  = cell(length(exp.participants),1); %pre-allocate
feat_AAC_Lv          = cell(length(exp.participants),1); %pre-allocate
feat_PACcan_Lv       = cell(length(exp.participants),1); %pre-allocate
feat_PACglm_Lv       = cell(length(exp.participants),1); %pre-allocate
feat_PACoz_Lv        = cell(length(exp.participants),1); %pre-allocate
feat_PACplv_Lv       = cell(length(exp.participants),1); %pre-allocate
feat_PACtort_Lv      = cell(length(exp.participants),1); %pre-allocate
feat_times_Lv        = cell(length(exp.participants),1); %pre-allocate

feat_lowFreq_amp_Riv = cell(length(exp.participants),1); %pre-allocate
feat_AAC_Riv         = cell(length(exp.participants),1); %pre-allocate
feat_PACcan_Riv      = cell(length(exp.participants),1); %pre-allocate
feat_PACglm_Riv      = cell(length(exp.participants),1); %pre-allocate
feat_PACoz_Riv       = cell(length(exp.participants),1); %pre-allocate
feat_PACplv_Riv      = cell(length(exp.participants),1); %pre-allocate
feat_PACtort_Riv     = cell(length(exp.participants),1); %pre-allocate
feat_times_Riv       = cell(length(exp.participants),1); %pre-allocate

feat_lowFreq_amp_Rv  = cell(length(exp.participants),1); %pre-allocate
feat_AAC_Rv          = cell(length(exp.participants),1); %pre-allocate
feat_PACcan_Rv       = cell(length(exp.participants),1); %pre-allocate
feat_PACglm_Rv       = cell(length(exp.participants),1); %pre-allocate
feat_PACoz_Rv        = cell(length(exp.participants),1); %pre-allocate
feat_PACplv_Rv       = cell(length(exp.participants),1); %pre-allocate
feat_PACtort_Rv      = cell(length(exp.participants),1); %pre-allocate
feat_times_Rv        = cell(length(exp.participants),1); %pre-allocate

for i_part = 1:length(exp.participants)

    % Select single subject data   
    partdata_Liv = data_out_Liv{i_part}; % chan x time x trial
    partdata_Lv = data_out_Lv{i_part}; % chan x time x trial
    partdata_Riv = data_out_Riv{i_part}; % chan x time x trial
    partdata_Rv = data_out_Rv{i_part}; % chan x time x trial

    % Get low frequency peaks defined by fooof
    partpk_Liv = alpha_trialpk_Liv{i_part}'; % chan x trial
    partpk_Lv = alpha_trialpk_Lv{i_part}'; % chan x trial
    partpk_Riv = alpha_trialpk_Riv{i_part}'; % chan x trial
    partpk_Rv = alpha_trialpk_Rv{i_part}'; % chan x trial


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute features for all the channels. %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through list of electrodes
    for ii = 1:length(select_elect)
        
        ichan = select_elect(ii); %for selecting electrodes

        % Get data from one brain electrode & transpose so trials x times
        singleChanEEG_Liv = squeeze(partdata_Liv(ichan,:,:))';
        singleChanEEG_Lv = squeeze(partdata_Lv(ichan,:,:))';
        singleChanEEG_Riv = squeeze(partdata_Riv(ichan,:,:))';
        singleChanEEG_Rv = squeeze(partdata_Rv(ichan,:,:))';
        
        % Add some padding before and after data
        singleChanEEG_Liv = ft_preproc_padding(singleChanEEG_Liv, 'zero', 600);
        singleChanEEG_Lv = ft_preproc_padding(singleChanEEG_Lv, 'zero', 600);
        singleChanEEG_Riv = ft_preproc_padding(singleChanEEG_Riv, 'zero', 600);
        singleChanEEG_Rv = ft_preproc_padding(singleChanEEG_Rv, 'zero', 600);
        
        
        % Get peak low freq for filtering
        lowFreq_Liv = squeeze(nanmedian(partpk_Liv(ii,:),2));
        lowFreq_Lv = squeeze(nanmedian(partpk_Lv(ii,:),2));
        lowFreq_Riv = squeeze(nanmedian(partpk_Riv(ii,:),2));
        lowFreq_Rv = squeeze(nanmedian(partpk_Rv(ii,:),2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract low-frequency phase & amplitude

        % Filter - Liv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Liv, Fs, [lowFreq_Liv-1 lowFreq_Liv+1],...
            filt_order, 'but', 'twopass', 'no');
        % Instant phase
        phase_filt_low_Liv(:,:) = ft_preproc_hilbert(filt, 'angle'); 
        % Instant amplitude
        amp_filt_low_Liv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Lv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Lv, Fs, [lowFreq_Lv-1 lowFreq_Lv+1],...
            filt_order, 'but', 'twopass', 'no');
        % Instant phase
        phase_filt_low_Lv(:,:) = ft_preproc_hilbert(filt, 'angle'); 
        % Instant amplitude
        amp_filt_low_Lv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Riv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Riv, Fs, [lowFreq_Riv-1 lowFreq_Riv+1],...
            filt_order, 'but', 'twopass', 'no');
        % Instant phase
        phase_filt_low_Riv(:,:) = ft_preproc_hilbert(filt, 'angle'); 
        % Instant amplitude
        amp_filt_low_Riv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Rv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Rv, Fs, [lowFreq_Rv-1 lowFreq_Rv+1],...
            filt_order, 'but', 'twopass', 'no');
        % Instant phase
        phase_filt_low_Rv(:,:) = ft_preproc_hilbert(filt, 'angle'); 
        % Instant amplitude
        amp_filt_low_Rv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract high-frequency phase & amplitude

        % Center amplitude frequency method
        Af1 = round(highFreq(1) -(mean(highFreq)/2.5));
        Af2 = round(highFreq(2) +(mean(highFreq)/2.5));
        
        % Max phase method
    %     Af1 = highFreq(1) - 1.5*max(lowFreq_Liv);
    %     Af2 = highFreq(2) + 1.5*max(lowFreq_Liv);

        % Filter - Liv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Liv, Fs, [Af1 Af2], filt_order,...
            'but', 'twopass', 'no');
        % Instant phase
        phase_filt_high_Liv(:,:) = ft_preproc_hilbert(filt, 'angle');
        % Instant amplitude
        amp_filt_high_Liv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Lv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Lv, Fs, [Af1 Af2], filt_order,...
            'but', 'twopass', 'no');
        % Instant phase
        phase_filt_high_Lv(:,:) = ft_preproc_hilbert(filt, 'angle');
        % Instant amplitude
        amp_filt_high_Lv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Riv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Riv, Fs, [Af1 Af2], filt_order,...
            'but', 'twopass', 'no');
        % Instant phase
        phase_filt_high_Riv(:,:) = ft_preproc_hilbert(filt, 'angle');
        % Instant amplitude
        amp_filt_high_Riv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        % Filter - Rv
        [filt] = ft_preproc_bandpassfilter(singleChanEEG_Rv, Fs, [Af1 Af2], filt_order,...
            'but', 'twopass', 'no');
        % Instant phase
        phase_filt_high_Rv(:,:) = ft_preproc_hilbert(filt, 'angle');
        % Instant amplitude
        amp_filt_high_Rv(:,:) = ft_preproc_hilbert(filt, 'abs');
        clear filt
        
        clear Af1 Af2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Set time window(s) of analysis
        % convert time window to logical
        timelim = times>=timewindow(1) & times<=timewindow(2);
        Ntime = length(times(timelim)); %length of analysis window

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FEATURE COMPUTATION - Liv
        ntrials = size(singleChanEEG_Liv,1); %number of trials in condition
        MI_tort_Liv    = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_oz_Liv      = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_can_Liv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_glm_Liv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_plv_Liv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        AAC_matrix_Liv = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        times_comp_Liv = NaN(length(1:(Ntime-winSize)/winStep),1); %pre-allocate
        
        % different # of trials per condition
        for trial = 1:ntrials
            phase_dataL = squeeze(phase_filt_low_Liv(trial,timelim));
            amp_dataH = squeeze(amp_filt_high_Liv(trial,timelim));
            amp_dataL = squeeze(amp_filt_low_Liv(trial,timelim));
            
            for ij = 1:(Ntime-winSize)/winStep
                
                % Get current time points for analysis
                tt = [(ij-1)*winStep+1:(ij-1)*winStep+winSize];
                
                % Mid-window for time
                times_comp_Liv(ij) = tt(round(winSize/2));
                
                %% PAC computation
%                 MI_tort_Liv(trial,ij) = calc_MI_tort(phase_dataL(tt),amp_dataH(tt),18);
                [MI_tort_Liv(trial,ij), ~] = eeg_klmi_OrientExo_ML(phase_dataL(tt),amp_dataH(tt),18);
                MI_oz_Liv(trial,ij) = calc_MI_ozkurt(phase_dataL(tt),amp_dataH(tt));
                MI_plv_Liv(trial,ij) = cohen_PLV(phase_dataL(tt),amp_dataH(tt));
                MI_can_Liv(trial,ij) = calc_MI_canolty(phase_dataL(tt),amp_dataH(tt)); 
                [MI_glm_Liv(trial,ij), ~] = eeg_glm_OrientExo_ML(phase_dataL(tt)',amp_dataH(tt)');
                
                %% AAC computation
                AAC = corr(amp_dataL(tt)',amp_dataH(tt)','type','s');              
                AAC_matrix_Liv(trial,ij) = AAC; % Add to matrix

                clear MI AAC tt
            end
            
            clear ij phase_dataL amp_dataH amp_dataL
        end
        clear trial ntrials
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FEATURE COMPUTATION - Lv
        ntrials = size(singleChanEEG_Lv,1); %number of trials in condition
        MI_tort_Lv    = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_oz_Lv      = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_can_Lv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_glm_Lv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_plv_Lv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        AAC_matrix_Lv = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        times_comp_Lv = NaN(length(1:(Ntime-winSize)/winStep),1); %pre-allocate
        
        % different # of trials per condition
        for trial = 1:ntrials
            
            phase_dataL = squeeze(phase_filt_low_Lv(trial,timelim));
            amp_dataH = squeeze(amp_filt_high_Lv(trial,timelim));
            amp_dataL = squeeze(amp_filt_low_Lv(trial,timelim));
            
            for ij = 1:(Ntime-winSize)/winStep
                
                % Get current time points for analysis
                tt = [(ij-1)*winStep+1:(ij-1)*winStep+winSize];
                
                % Mid-window for time
                times_comp_Lv(ij) = tt(round(winSize/2));
                
                %% PAC computation
                [MI_tort_Lv(trial,ij), ~] = eeg_klmi_OrientExo_ML(phase_dataL(tt),amp_dataH(tt),18);
                MI_oz_Lv(trial,ij) = calc_MI_ozkurt(phase_dataL(tt),amp_dataH(tt));
                MI_plv_Lv(trial,ij) = cohen_PLV(phase_dataL(tt),amp_dataH(tt));
                MI_can_Lv(trial,ij) = calc_MI_canolty(phase_dataL(tt),amp_dataH(tt)); 
                [MI_glm_Lv(trial,ij), ~] = eeg_glm_OrientExo_ML(phase_dataL(tt)',amp_dataH(tt)');
                
                %% AAC computation
                AAC = corr(amp_dataL(tt)',amp_dataH(tt)','type','s');              
                AAC_matrix_Lv(trial,ij) = AAC; % Add to matrix

                clear MI AAC tt
            end
            
            clear ij phase_dataL amp_dataH amp_dataL
        end
        clear trial ntrials
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FEATURE COMPUTATION - Riv
        ntrials = size(singleChanEEG_Riv,1); %number of trials in condition
        MI_tort_Riv    = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_oz_Riv      = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_can_Riv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_glm_Riv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_plv_Riv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        AAC_matrix_Riv = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        times_comp_Riv = NaN(length(1:(Ntime-winSize)/winStep),1); %pre-allocate
        
        % different # of trials per condition
        for trial = 1:ntrials
            phase_dataL = squeeze(phase_filt_low_Riv(trial,timelim));
            amp_dataH = squeeze(amp_filt_high_Riv(trial,timelim));
            amp_dataL = squeeze(amp_filt_low_Riv(trial,timelim));
            
            for ij = 1:(Ntime-winSize)/winStep
                
                % Get current time points for analysis
                tt = [(ij-1)*winStep+1:(ij-1)*winStep+winSize];
                
                % Mid-window for time
                times_comp_Riv(ij) = tt(round(winSize/2));
                
                %% PAC computation
                [MI_tort_Riv(trial,ij), ~] = eeg_klmi_OrientExo_ML(phase_dataL(tt),amp_dataH(tt),18);
                MI_oz_Riv(trial,ij) = calc_MI_ozkurt(phase_dataL(tt),amp_dataH(tt));
                MI_plv_Riv(trial,ij) = cohen_PLV(phase_dataL(tt),amp_dataH(tt));
                MI_can_Riv(trial,ij) = calc_MI_canolty(phase_dataL(tt),amp_dataH(tt)); 
                [MI_glm_Riv(trial,ij), ~] = eeg_glm_OrientExo_ML(phase_dataL(tt)',amp_dataH(tt)');
                
                %% AAC computation
                AAC = corr(amp_dataL(tt)',amp_dataH(tt)','type','s');              
                AAC_matrix_Riv(trial,ij) = AAC; % Add to matrix

                clear MI AAC tt
            end
            
            clear ij phase_dataL amp_dataH amp_dataL
        end
        clear trial ntrials
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FEATURE COMPUTATION - Rv
        ntrials = size(singleChanEEG_Rv,1); %number of trials in condition
        MI_tort_Rv    = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_oz_Rv      = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_can_Rv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_glm_Rv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        MI_plv_Rv     = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        AAC_matrix_Rv = NaN(ntrials,length(1:(Ntime-winSize)/winStep)); %pre-allocate
        times_comp_Rv = NaN(length(1:(Ntime-winSize)/winStep),1); %pre-allocate
        
        % different # of trials per condition
        for trial = 1:ntrials
            
            phase_dataL = squeeze(phase_filt_low_Rv(trial,timelim));
            amp_dataH = squeeze(amp_filt_high_Rv(trial,timelim));
            amp_dataL = squeeze(amp_filt_low_Rv(trial,timelim));
            
            for ij = 1:(Ntime-winSize)/winStep
                
                % Get current time points for analysis
                tt = [(ij-1)*winStep+1:(ij-1)*winStep+winSize];
                
                % Mid-window for time
                times_comp_Rv(ij) = tt(round(winSize/2));
                
                %% PAC computation
                [MI_tort_Rv(trial,ij), ~] = eeg_klmi_OrientExo_ML(phase_dataL(tt),amp_dataH(tt),18);
                MI_oz_Rv(trial,ij) = calc_MI_ozkurt(phase_dataL(tt),amp_dataH(tt));
                MI_plv_Rv(trial,ij) = cohen_PLV(phase_dataL(tt),amp_dataH(tt));
                MI_can_Rv(trial,ij) = calc_MI_canolty(phase_dataL(tt),amp_dataH(tt)); 
                [MI_glm_Rv(trial,ij), ~] = eeg_glm_OrientExo_ML(phase_dataL(tt)',amp_dataH(tt)');
                
                %% AAC computation
                AAC = corr(amp_dataL(tt)',amp_dataH(tt)','type','s');              
                AAC_matrix_Rv(trial,ij) = AAC; % Add to matrix

                clear MI AAC tt
            end
            
            clear ij phase_dataL amp_dataH amp_dataL
        end
        clear trial ntrials

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save features from each electrode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        feat_times_Liv{i_part,1} = times_comp_Liv;
        feat_times_Lv{i_part,1} = times_comp_Lv;
        feat_times_Riv{i_part,1} = times_comp_Riv;
        feat_times_Rv{i_part,1} = times_comp_Rv;
        
        feat_lowFreq_amp_Liv{i_part,1}(:,ii,:) = squeeze(amp_filt_low_Liv(:,timelim))';
        feat_lowFreq_amp_Lv{i_part,1}(:,ii,:) = squeeze(amp_filt_low_Lv(:,timelim))';
        feat_lowFreq_amp_Riv{i_part,1}(:,ii,:) = squeeze(amp_filt_low_Riv(:,timelim))';
        feat_lowFreq_amp_Rv{i_part,1}(:,ii,:) = squeeze(amp_filt_low_Rv(:,timelim))';
        
        feat_AAC_Liv{i_part,1}(:,ii,:) = AAC_matrix_Liv;
        feat_AAC_Lv{i_part,1}(:,ii,:) = AAC_matrix_Lv;
        feat_AAC_Riv{i_part,1}(:,ii,:) = AAC_matrix_Riv;
        feat_AAC_Rv{i_part,1}(:,ii,:) = AAC_matrix_Rv;
        
        feat_PACcan_Liv{i_part,1}(:,ii,:) = MI_can_Liv;
        feat_PACcan_Lv{i_part,1}(:,ii,:) = MI_can_Lv;
        feat_PACcan_Riv{i_part,1}(:,ii,:) = MI_can_Riv;
        feat_PACcan_Rv{i_part,1}(:,ii,:) = MI_can_Rv;
        
        feat_PACglm_Liv{i_part,1}(:,ii,:) = MI_glm_Liv;
        feat_PACglm_Lv{i_part,1}(:,ii,:) = MI_glm_Lv;
        feat_PACglm_Riv{i_part,1}(:,ii,:) = MI_glm_Riv;
        feat_PACglm_Rv{i_part,1}(:,ii,:) = MI_glm_Rv;
        
        feat_PACoz_Liv{i_part,1}(:,ii,:) = MI_oz_Liv;
        feat_PACoz_Lv{i_part,1}(:,ii,:) = MI_oz_Lv;
        feat_PACoz_Riv{i_part,1}(:,ii,:) = MI_oz_Riv;
        feat_PACoz_Rv{i_part,1}(:,ii,:) = MI_oz_Rv;
        
        feat_PACplv_Liv{i_part,1}(:,ii,:) = MI_plv_Liv;
        feat_PACplv_Lv{i_part,1}(:,ii,:) = MI_plv_Lv;
        feat_PACplv_Riv{i_part,1}(:,ii,:) = MI_plv_Riv;
        feat_PACplv_Rv{i_part,1}(:,ii,:) = MI_plv_Rv;
        
        feat_PACtort_Liv{i_part,1}(:,ii,:) = MI_tort_Liv;
        feat_PACtort_Lv{i_part,1}(:,ii,:) = MI_tort_Lv;
        feat_PACtort_Riv{i_part,1}(:,ii,:) = MI_tort_Riv;
        feat_PACtort_Rv{i_part,1}(:,ii,:) = MI_tort_Rv;
        

        clear AAC_matrix* MI_matrix* singleChanEEG* phase_filt_low*...
            amp_filt_low* phase_filt_high* amp_filt_high* MI* times_comp*...
            lowFreq*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %channel loop
    
    clear ii ichan partdata* partpk*

end
clear i_part beta


% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////

%% Save with version for large files
times = times_out; %before padding

save([pwd '/feats_v2_times' num2str(length(feat_times_Rv{1})) '_' exp.settings '.mat'],...
    'timewindow','times','select_elect','filt_order','winSize','winStep',...
    '-regexp','feat*','-v7.3');

% save([pwd '/feats_times' num2str(length(feat_times_Rv{1})) '_' exp.settings '.mat'],...
%     '-regexp','feat*','-append','-v7.3');


% /////////////////////////////////////////////////////////////////////////
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////

clear Fs filt_order highFreq winSize winStep method Ntime


feat_plot = squeeze(feat_AAC_Liv{1}(:,1,:));
time_plot = times(timelim);

figure
plot(time_plot(feat_times_Liv{1,1}),feat_plot)


featavg_plot = squeeze(mean(feat_AAC_Riv{1}(:,1,:),1));

figure
plot(time_plot(feat_times_Liv{1,1}),featavg_plot)




