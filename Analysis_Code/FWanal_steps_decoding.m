%% Housekeeping

% Clears the workspace and closes all figure windows
clear variables;
% close all;

% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% /////////////////////////////////////////////////////////////////////////
%% General Settings

% Enter the name of the study (for labeling saved decoding results files)
% study_name = 'AAC';
% study_name = 'lowFreq_amp';
study_name = 'PACoz';
% study_name = 'PACplv';
% study_name = 'PACtort';
% study_name = 'PACglm';

% Base directory path (where single subject EEG datasets and channel locations files are stored)
maindir = [pwd '\decode_v2\' exp.settings '\']; %gets written over below
bdir = [maindir study_name '\'];

% Output directory (where decoding results will be saved)
output_dir = [bdir '\Decode_Results\']; %for classification 

% Which group-level statistical analysis method was used??
% 1 = Global null and population prevalence tests based on the minimum statistic
% 2 = Global null testing using paired-samples t tests
group_level_analysis_method = 2; 

% Which decoding result? 1 = right targets, 2 = left targets
file_n = 1;

% -------------------------------------------------------------------------
%% Appropriate labels

if strcmpi(study_name,'lowFreq_amp')
    % lowFreq_amp data
    filename{1} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVM_LIBSVM_DCGR_NI vs. R_I.mat'];
    filename{2} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVM_LIBSVM_DCGL_NI vs. L_I.mat'];
else  
    % All other feature
    filename{1} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVM_LIBSVM_DCGR_NI vs. R_I.mat'];
    filename{2} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVM_LIBSVM_DCGL_NI vs. L_I.mat'];
end
    

% -------------------------------------------------------------------------
%% Location to save plots

saveFig = [bdir 'Figures\GROUPRES_SVM\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

% -------------------------------------------------------------------------
%% Load file
group_results_file = [output_dir filename{file_n}];

% Load the data file to get ANALYSIS parameters
load(group_results_file);

clear group_results_file


% transfer info
FW.chanlocs = FW_ANALYSIS.chanlocs;
FW.chaninfo = FW_ANALYSIS.chaninfo;

% -------------------------------------------------------------------------
%% Change time steps to x-axis scale - help choose steps from figure

step_scale = ANALYSIS.xaxis_scale;

if strcmpi(study_name,'lowFreq_amp')
    step_scale(2,1:end) = 0 : ANALYSIS.step_width_ms : 999;
else
    step_scale(2,1:end) = 125 : ANALYSIS.step_width_ms : 874;
end

% Change labels for plots and save names
step_scale_label = ANALYSIS.xaxis_scale;

if strcmpi(study_name,'lowFreq_amp')
    step_scale_label(2,1:end) = (0 : ANALYSIS.step_width_ms : 999) - ANALYSIS.pointzero;
else
    step_scale_label(2,1:end) = (125 : ANALYSIS.step_width_ms : 874) - ANALYSIS.pointzero;
end

% Enter the time-steps for statistical testing of feature weights (e.g. [4 6 7 10])
%my version is in FW to not overwrite the original
% Example: fw.steps_for_testing = [5:10]
% FW.steps_for_testing = [1:50]; %all of the steps

% FW.steps_for_testing = 6:31; %low_freq - L
% FW.steps_for_testing = 38:50; %low_freq - L
% FW.steps_for_testing = 2:5; %low_freq - R
% FW.steps_for_testing = 6:17; %low_freq - R
% FW.steps_for_testing = 27:48; %low_freq - R

% FW.steps_for_testing = 1; %PACoz - L
% FW.steps_for_testing = 7:8; %PACoz - R
FW.steps_for_testing = 34; %PACoz - R
FW.disp_steps = FW.steps_for_testing;


% Add steps to plot names
plotname{1} = [study_name ': Non-Info vs. Info, R Target (' ...
    num2str(step_scale_label(2,FW.steps_for_testing(1))) '-' num2str(step_scale_label(2,FW.steps_for_testing(end))) ' ms)'];
plotname{2} = [study_name ': Non-Info vs. Info, L Target (' ...
    num2str(step_scale_label(2,FW.steps_for_testing(1))) '-' num2str(step_scale_label(2,FW.steps_for_testing(end))) ' ms)'];

% Names to save plots
figname{1} = [study_name '_' num2str(step_scale_label(2,FW.steps_for_testing(1))) '-' ...
    num2str(step_scale_label(2,FW.steps_for_testing(end))) 'ms_FW_SVM_R_NIvsI.fig'];
figname{2} = [study_name '_' num2str(step_scale_label(2,FW.steps_for_testing(1))) '-' ...
    num2str(step_scale_label(2,FW.steps_for_testing(end))) 'ms_FW_SVM_L_NIvsI.fig'];

% -------------------------------------------------------------------------
%% Get data from each participant

for sbj = 1:ANALYSIS.nsbj
           
    % results are stored in:
    % RESULTS.feature_weights{analysis}{time_step, cross_val_step, repetition_step}(feature_number, feature_weight, absolute_feature_weight)
    %
    % Feature weights are then averaged. Attention - averaged absolute
    % feature weights are not just | average feature weights | (!!!)
        
    for steps = 1:size(ANALYSIS.RES.feature_weights{sbj}, 1)

        % Preallocate summed feature weights matrices
        sz = size(ANALYSIS.RES.feature_weights{sbj}{steps,1,1}, 1);
        sum_step_fw(sz, 1) = zeros;
        sum_step_abs(sz, 1) = zeros;

        for cross_val = 1:size(ANALYSIS.RES.feature_weights{sbj}, 2)

            for rep = 1:size(ANALYSIS.RES.feature_weights{sbj}, 3)

                if ANALYSIS.fw.corrected == 0 % if analysing uncorrected feature weights

                    temp_step_fw = ANALYSIS.RES.feature_weights{1, sbj}{steps, cross_val, rep}(:,2);
                    temp_step_abs = ANALYSIS.RES.feature_weights{1, sbj}{steps, cross_val, rep}(:,3); 

                elseif ANALYSIS.fw.corrected == 1 % if analysing corrected feature weights

                    temp_step_fw = ANALYSIS.RES.feature_weights_corrected{1, sbj}{steps, cross_val, rep}(:,2);
                    temp_step_abs = ANALYSIS.RES.feature_weights_corrected{1, sbj}{steps, cross_val, rep}(:,3); 

                end % of if ANALYSIS.fw.corrected

                FW.ALL_FW{sbj, steps, cross_val, rep} = temp_step_fw; % feature weights
                FW.ALL_FW{sbj, steps, cross_val, rep}(:,2) = temp_step_abs; % absolute feature weights

                sum_step_fw = sum_step_fw + temp_step_fw;
                clear temp_step_fw;
                sum_step_abs = sum_step_abs + temp_step_abs;
                clear temp_step_abs;

            end % of for rep

        end % of for cross_val

        % Average feature weights
        FW.ALL_FW_STEP{sbj,steps}(:,1) = (sum_step_fw ./ (size(ANALYSIS.RES.feature_weights{sbj}, 2) + size(ANALYSIS.RES.feature_weights{sbj}, 3)));
        FW.ALL_FW_STEP{sbj,steps}(:,2) = (sum_step_abs ./ (size(ANALYSIS.RES.feature_weights{sbj}, 2) + size(ANALYSIS.RES.feature_weights{sbj}, 3)));

        clear sum_step_abs;
        clear sum_step_fw;

    end % of for steps

    %% Average FWs Across Time WITHIN an Analysis Time Window (For Temporal/Spatio-Temporal Decoding)
    % not required for spatial decoding, because that is based on average
    % of time step
    
    if ANALYSIS.stmode == 3 % Spatio-temporal decoding
                   
       for steps = 1:size(FW.ALL_FW_STEP, 2)
                    
            start_point = 1;
            
            for withinstep = 1:ANALYSIS.nchannels
                    
                temp_fw = FW.ALL_FW_STEP{sbj, steps}(start_point:(start_point + ANALYSIS.window_width - 1), 1);
                temp_abs = FW.ALL_FW_STEP{sbj, steps}(start_point:(start_point + ANALYSIS.window_width - 1), 2);
                       
                m_temp_fw = mean(temp_fw);
                m_temp_abs = mean(temp_abs);
                       
                FW.ALL_FW_STEP_CHANNEL{sbj, steps}(withinstep, 1) = m_temp_fw;
                FW.ALL_FW_STEP_CHANNEL{sbj, steps}(withinstep, 2) = m_temp_abs;
                       
                start_point = start_point + ANALYSIS.window_width;
                        
                clear temp_fw; 
                clear temp_abs;
                clear m_temp_fw; 
                clear m_temp_abs;
                                                
            end % of for withinstep
                    
        end % of for steps
    
    else
        
        FW.ALL_FW_STEP_CHANNEL = FW.ALL_FW_STEP;
        
    end % of if ANALYSIS.stmode
    
end % of for sbj


%% Transform Absolute Feature Weights For Each Channel Into Z-Scores
% this is done separately for each analysis (as comparing feature weights
% between analyses is not meaningful, because they are relative to others for the same analysis)

for sbj = 1:size(FW.ALL_FW_STEP_CHANNEL, 1)

    for steps = 1:size(FW.ALL_FW_STEP_CHANNEL, 2)
           
        temp(:,1) = FW.ALL_FW_STEP_CHANNEL{sbj, steps}(:,2);
        fw_all_abs(sbj, steps, :) = temp;
        
        % convert to z-scores
        temp_z = zscore(temp);
        fw_all_z(sbj, steps, :) = temp_z;
        
        % third column contains z-scores of absolute feature weights
        FW.ALL_FW_STEP_CHANNEL{sbj, steps}(:,3) = temp_z;
        
        clear temp;
        clear temp_z;
                
    end % of for steps
end % of for sbj

% Copy to ANALYSIS structure
FW.ALL_Z = fw_all_z;
clear fw_all_z;
FW.ALL_ABS = fw_all_abs;
clear fw_all_abs;


%% Average Z-Scores Across Participants For Each Step
% Generates FW_ANALYSIS.AVERAGE_Z(step, channel)

mean_fw_all_z(:,:) = mean(FW.ALL_Z, 1);
FW.AVERAGE_Z = mean_fw_all_z;

% Analysis time windows for matrix display
FW.AVERAGE_Z_DISP = FW.AVERAGE_Z(FW.disp_steps(1):FW.disp_steps(end), :);

% Analysis time windows for heat map displays
for stat_step = 1:size(FW.steps_for_testing, 2)
    
    FW.AVERAGE_Z_HEATS(stat_step, :) = FW.AVERAGE_Z(FW.disp_steps(stat_step), :);
    
end % of for stat_step


%% Average Relevant Time Steps For Stats/Display

sum_fw = [];

for steps = 1:size(FW.steps_for_testing, 2)
    
    % Extract absolute and Z-scored feature weights at each analysis step
    temp_fw_select(:,:) = FW.ALL_ABS(:, FW.steps_for_testing(steps), :);
    temp_fw_select_z(:,:) = FW.ALL_Z(:, FW.steps_for_testing(steps), :);
    
    % Preallocate the feature weight sum matrices for the first step
    if steps == 1
        
        sz = size(temp_fw_select);
        sum_fw(sz(1), sz(2)) = zeros;
        sum_fw_z(sz(1), sz(2)) = zeros;
        
    end % of if steps
    
    sum_fw = sum_fw + temp_fw_select;
    sum_fw_z = sum_fw_z + temp_fw_select_z;
    
    clear temp_fw_pre;
    clear temp_fw_pre_z;
    
end % of for steps

% Average absolute and Z-scored feature weights
FW.AVERAGESTEPS_SELECT_FW_ABS = sum_fw ./ size(FW.steps_for_testing,2);
FW.AVERAGESTEPS_SELECT_FW_Z = sum_fw_z ./ size(FW.steps_for_testing,2);

FW.AVERAGESTEPS_SELECT_FW_ABS_MEAN = mean(FW.AVERAGESTEPS_SELECT_FW_ABS,1);
FW.AVERAGESTEPS_SELECT_FW_Z_MEAN = mean(FW.AVERAGESTEPS_SELECT_FW_Z,1);



%% T Tests For FWs - Averaged Selected Time Bins & Single Selected Time Bins
% Output is one matrix for each analysis, once with correction for multiple
% comparisons ( = features/channels) and once without (p < critical value)

% T-test for all single analysis time steps
for p_corr = 1:2 % run for corrected/uncorrected
    
    if p_corr == 1 % Uncorrected for multiple comparisons
        
        p_crit = ANALYSIS.pstats; % Uncorrected alpha level

        for steps = 1:size(FW.steps_for_testing, 2)
            
            for channel = 1:size(FW.ALL_Z, 3)
                
                temp_z = FW.ALL_Z(:, FW.steps_for_testing(steps), channel);

                % Run a one-sample t-test on Z-scored feature weights
                
                if ANALYSIS.fw.use_robust == 0 % Student's t test
                    
                    [h,p] = ttest(temp_z, 0, p_crit, ANALYSIS.fw.ttest_tail); 
                
                elseif ANALYSIS.fw.use_robust == 1 % Yuen's t
                    
                    zero_data_temp = zeros(length(temp_z), 1); % Make vector of zeroes for single-sample comparison
                    [h,p, ~, ~, ~, ~, ~, ~] = yuend_ttest(temp_z, zero_data_temp, ...
                        'percent', ANALYSIS.fw.trimming, 'alpha', ANALYSIS.pstats, ...
                        'tail', ANALYSIS.fw.ttest_tail);

                end % of if ANALYSIS.fw.use_robust
                
                h_matrix_z(channel, 1) = h; % Stores significant/non-significant decisions
                p_matrix_z(channel, 1) = p; % Stores p-values

                clear temp_z;
                clear h;
                clear p;

            end % of for channel

            % Copy t-test results into FW_ANALYSIS structure
            FW.p_matrix_z_uncorr{steps} = p_matrix_z; 
            FW.p_matrix_z_uncorr_label = FW.steps_for_testing;
            FW.h_matrix_z_uncorr{steps} = h_matrix_z;
            clear p_matrix_z;
            clear h_matrix_z;
        
        end % of for steps

        elseif p_corr == 2 % Corrected for multiple comparisons
            
            switch ANALYSIS.fw.multcompstats % Selects a multiple comparisons correction 
                
                case 1 % Bonferroni correction
                    
                    fprintf('\nPerforming corrections for multiple comparisons (Bonferroni)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; % Copy same p-values as uncorrected
                        [FW_MCC] = multcomp_bonferroni(FW.p_matrix_z_uncorr{steps}, 'alpha', ANALYSIS.fw.pstats);
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        
                    end % of for steps
            
                case 2 % Holm-Bonferroni correction
                    
                    fprintf('\nPerforming corrections for multiple comparisons (Holm-Bonferroni)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; % Copy same p-values as uncorrected
                        [FW_MCC] = multcomp_holm_bonferroni(FW.p_matrix_z_uncorr{steps}, 'alpha', ANALYSIS.fw.pstats);
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        
                    end % of for steps
                    
                case 3 % strong FWER control permutation test
                    
                    fprintf('\nPerforming corrections for multiple comparisons (maximum statistic permutation test)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        real_decoding_fw = FW.ALL_Z(:,FW.steps_for_testing(steps),:); % Results matrix (subjects x channels)
                        real_decoding_fw = squeeze(real_decoding_fw); % Remove extra dimension defined by step number
                        fw_chance_level = zeros(size(real_decoding_fw, 1), size(real_decoding_fw, 2)); % Matrix of chance level values (zeros)
                        
                        [FW_MCC] = multcomp_blair_karniski_permtest(real_decoding_fw, fw_chance_level, ...
                            'alpha', ANALYSIS.fw.pstats, 'iterations', ANALYSIS.fw.n_iterations, ...
                            'use_yuen', ANALYSIS.fw.use_robust, 'percent', ANALYSIS.fw.trimming, ...
                            'tail', ANALYSIS.fw.ttest_tail);
                        
                        % Copy results to FW_ANALYSIS
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        FW.p_matrix_z_corr{steps} = FW_MCC.corrected_p;
                        
                    end % of for steps
                    
                    clear real_decoding_fw;
                    clear fw_chance_level;
                    
                case 4 % cluster-based permutation test (To be implemented in future, once we can easily set up neighborhood matrices)
                    
                    % Stopgap until cluster-based permutation testing is
                    % implemented (needs neighborhood matrix of electrode
                    % posititons).
                    fprintf('\nCluster-based multiple comparisons correction method not currently available... No correction was performed\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label;
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; 
                        FW.h_matrix_z_corr{steps} = FW.h_matrix_z_uncorr{steps};
                        
                    end % of for steps
                    
                case 5 % Generalised FWER control procedure
                    
                    fprintf('\nPerforming corrections for multiple comparisons (KTMS generalised FWER control)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        real_decoding_fw = FW.ALL_Z(:,FW.steps_for_testing(steps),:); % Results matrix (subjects x channels)
                        real_decoding_fw = squeeze(real_decoding_fw); % Remove extra dimension defined by step number
                        fw_chance_level = zeros(size(real_decoding_fw, 1), size(real_decoding_fw, 2)); % Matrix of chance level values (zeros)
                        
                        [FW_MCC] = multcomp_ktms(real_decoding_fw, fw_chance_level, ...
                            'alpha', ANALYSIS.fw.pstats, 'iterations', ANALYSIS.fw.n_iterations, ...
                            'ktms_u', ANALYSIS.fw.ktms_u, 'use_yuen', ANALYSIS.fw.use_robust, ...
                            'percent', ANALYSIS.fw.trimming, 'tail', ANALYSIS.fw.ttest_tail);
                        
                        % Copy results to FW_ANALYSIS
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        FW.p_matrix_z_corr{steps} = FW_MCC.corrected_p;
                        
                    end % of for steps
                    
                    clear real_decoding_fw;
                    clear fw_chance_level;
                    
                case 6 % Benjamini-Hochberg false discovery rate control
                    
                    fprintf('\nPerforming corrections for multiple comparisons (Benjamini-Hochberg FDR control)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; % Copy same p-values as uncorrected
                        [FW_MCC] = multcomp_fdr_bh(FW.p_matrix_z_uncorr{steps}, 'alpha', ANALYSIS.fw.pstats); 
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        
                    end % of for steps
                    
                case 7 % Benjamini-Krieger-Yekutieli false discovery rate control
                    
                    fprintf('\nPerforming corrections for multiple comparisons (Benjamini-Krieger-Yekutieli FDR control)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; % Copy same p-values as uncorrected
                        FW_MCC = multcomp_fdr_bky(FW.p_matrix_z_uncorr{steps}, 'alpha', ANALYSIS.fw.pstats);
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        
                    end % of for steps
                      
                case 8 % Benjamini-Yekutieli false discovery rate control
                    
                    fprintf('\nPerforming corrections for multiple comparisons (Benjamini-Yekutieli FDR control)\n\n');
                    
                    FW.p_matrix_z_corr_label = FW.p_matrix_z_uncorr_label; % Copy same labels as uncorrected
                    
                    for steps = 1:size(FW.steps_for_testing, 2)
                        
                        FW.p_matrix_z_corr{steps} = FW.p_matrix_z_uncorr{steps}; % Copy same p-values as uncorrected
                        FW_MCC = multcomp_fdr_by(FW.p_matrix_z_uncorr{steps}, 'alpha', ANALYSIS.fw.pstats);
                        FW.h_matrix_z_corr{steps} = FW_MCC.corrected_h;
                        
                    end % of for steps

            end % of ANALYSIS.fw.multcompstats switch
    end % of if p_corr
end % of for p_corr



%% T-Tests For Averaged Analysis Time Window

for p_corr = 1:2 % run for corrected/uncorrected
    
    if p_corr == 1 % Uncorrected for multiple comparisons
        
        p_crit = ANALYSIS.pstats; % Uncorrected alpha level

        for channel = 1:size(FW.ALL_Z, 3)

            temp = FW.AVERAGESTEPS_SELECT_FW_Z(:, channel);
%             temp = FW_ANALYSIS.AVERAGESTEPS_SELECT_FW_Z(channel, :);

            if ANALYSIS.fw.use_robust == 0 % Student's t test
                
                [h,p] = ttest(temp, 0, p_crit, ANALYSIS.fw.ttest_tail); 
            
            elseif ANALYSIS.fw.use_robust == 1 % Yuen's t test
            
                zero_data_temp = zeros(length(temp), 1); % Make vector of zeroes for single-sample comparison
                [h,p, ~, ~, ~, ~, ~, ~] = yuend_ttest(temp, zero_data_temp, 'percent', ANALYSIS.fw.trimming, 'alpha', ANALYSIS.pstats, 'tail', ANALYSIS.fw.ttest_tail);
                
            end % of if ANALYSIS.fw.use_robust
            
            h_matrix_z(channel,1) = h; % Stores significant/non-significant decisions
            p_matrix_z(channel,1) = p; % Stores p-values

        end  % of for channel

        % Copy t-test results into FW_ANALYSIS structure
        FW.p_matrix_z_averagestep_uncorr = p_matrix_z; 
        FW.p_matrix_z_averagestep_uncorr_label = FW.steps_for_testing;
        FW.h_matrix_z_averagestep_uncorr = h_matrix_z;
        clear p_matrix_z;
        clear h_matrix_z;
        
    elseif p_corr == 2 % Corrected for multiple comparisons
        
        % Multiple comparisons corrections
        switch ANALYSIS.fw.multcompstats % Selects a multiple comparisons correction 
            
                case 1 % Bonferroni correction          
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; % Copy from uncorrected ver.
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label; % Copy from uncorrected ver.
                    [FW_MCC] = multcomp_bonferroni(FW.p_matrix_z_averagestep_uncorr, 'alpha', ANALYSIS.fw.pstats); 
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                     
                case 2 % Holm-Bonferroni correction
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; % Copy from uncorrected ver.
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label; % Copy from uncorrected ver.
                    [FW_MCC] = multcomp_holm_bonferroni(FW.p_matrix_z_averagestep_uncorr, 'alpha', ANALYSIS.fw.pstats); 
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                    
                case 3 % Strong FWER control permutation test
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; % Copy from uncorrected ver.
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label;
                    real_decoding_fw = FW.AVERAGESTEPS_SELECT_FW_Z(:,:); % Results matrix (subjects x channels)
                    fw_chance_level = zeros(size(real_decoding_fw, 1), size(real_decoding_fw, 2), size(real_decoding_fw, 3)); % Matrix of chance level values (zeros)
                    
                    [FW_MCC] = multcomp_blair_karniski_permtest(real_decoding_fw, fw_chance_level, ...
                        'alpha', ANALYSIS.fw.pstats, 'iterations', ANALYSIS.fw.n_iterations, ...
                        'use_yuen', ANALYSIS.fw.use_robust, 'percent', ANALYSIS.fw.trimming, ...
                        'tail', ANALYSIS.fw.ttest_tail);
                    
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                    FW.p_matrix_z_averagestep_corr = FW_MCC.corrected_p;
                    FW.p_matrix_z_averagestep_corr_label = FW_MCC.corrected_p;
                    clear real_decoding_fw;
                    clear fw_chance_level;
                    
                case 4 % Cluster-based permutation test (not yet available)
                    
                    fprintf('\nCluster-based correction not currently available... No correction was performed\n');
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; % Copy from uncorrected ver.
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label; % Copy from uncorrected ver.
                    FW.h_matrix_z_averagestep_corr = FW.h_matrix_z_averagestep_uncorr;
                                        
                case 5 % Generalised FWER control procedure
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; % Copy from uncorrected ver.
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label;
                    real_decoding_fw = FW.AVERAGESTEPS_SELECT_FW_Z(:,:); % Results matrix (subjects x channels)
                    fw_chance_level = zeros(size(real_decoding_fw, 1), size(real_decoding_fw, 2)); % Matrix of chance level values (zeros)
                    
                    [FW_MCC] = multcomp_ktms(real_decoding_fw, fw_chance_level, ...
                        'alpha', ANALYSIS.fw.pstats, 'iterations', ANALYSIS.fw.n_iterations, ...
                        'ktms_u', ANALYSIS.fw.ktms_u, 'use_yuen', ANALYSIS.fw.use_robust, ...
                        'percent', ANALYSIS.fw.trimming, 'tail', ANALYSIS.fw.ttest_tail);
                    
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                    clear real_decoding_fw;
                    clear fw_chance_level;

                case 6 % Benjamini-Hochberg false discovery rate control
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; 
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label;
                    [FW_MCC] = multcomp_fdr_bh(FW.p_matrix_z_averagestep_uncorr, 'alpha', ANALYSIS.fw.pstats);
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                    
                case 7 % Benjamini-Krieger-Yekutieli false discovery rate control
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; 
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label;
                    [FW_MCC] = multcomp_fdr_bky(FW.p_matrix_z_averagestep_uncorr, 'alpha', ANALYSIS.fw.pstats);
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;
                      
                case 8 % Benjamini-Yekutieli false discovery rate control
                    
                    FW.p_matrix_z_averagestep_corr = FW.p_matrix_z_averagestep_uncorr; 
                    FW.p_matrix_z_averagestep_corr_label = FW.p_matrix_z_averagestep_uncorr_label;
                    [FW_MCC] = multcomp_fdr_by(FW.p_matrix_z_averagestep_uncorr, 'alpha', ANALYSIS.fw.pstats);
                    FW.h_matrix_z_averagestep_corr = FW_MCC.corrected_h;

        end % of ANALYSIS.fw.multcompstats switch
    end % of if p_corr
end % of for p_corr


% ---------------------------------------------------------------------
%% Plots!!!

% my code for plotting
FWanal_decoding_plot(ANALYSIS,FW,saveFig,figname,plotname,study_name,file_n);




