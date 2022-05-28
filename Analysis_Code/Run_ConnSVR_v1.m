%% Housekeeping

% Clears the workspace and closes all figure windows
clear variables;
close all;

%% Add necessary folders to path list
addpath(genpath('C:\Users\ssshe\Documents\MathLab\Personal_Folders\Sarah\Matlab_Functions\DDTBOX'))
%stupid computer is getting fixed
% addpath(genpath('D:\ssshe\Documents\MathLab\Personal_Folders\Sarah\Matlab_Functions\DDTBOX'))

% /////////////////////////////////////////////////////////////////////////
%% Load processing settingsit wou
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data
% load([pwd '/BEH_' exp.settings '.mat']);
% 
% % load specific EEG dataset to make EEGLab happy
% EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Select Subject Datasets and Discrimination Groups (dcgs)

% Set the subject datasets on which to perform MVPA
% (splitting up because takes so long)
% sbj_todo = [1:7];
sbj_todo = [8:14];
% sbj_todo = [15:21];
% sbj_todo = [22:28];
% sbj_todo = [1:length(exp.participants)]; %not based on actual id #

% Enter the discrimination groups (dcgs) for decoding analyses. 
% Each discrimination group should be in a separate cell entry.
% Decoding analyses will be run for all dcgs listed here.
% e.g. dcgs_for_analyses{1} = [1];
% e.g. dcgs_for_analyses{2} = [3];
% Two discrimination groups can be entered when performing cross-condition decoding.
% (SVM trained using the first entry/dcg, tested on the second entry/dcg)
% e.g. dcgs_for_analyses{1} = [1, 2];

% 1 = L_NI    2 = L_I    3 = R_NI    4 = R_I

dcgs_for_analyses{1} = [1];
dcgs_for_analyses{2} = [2];
dcgs_for_analyses{3} = [3];
dcgs_for_analyses{4} = [4];


% Perform cross-condition decoding? 
% 0 = No / 1 = Yes
cross = 0;



%% Filepaths and Locations of Subject Datasets

% Enter the name of the study (for labeling saved decoding results files)
% study_name = 'corr';
% study_name = 'dwPLI';
study_name = 'iCOH';
% % study_name = 'PLI';
% study_name = 'PLV';
% % study_name = 'wPLI';

% Base directory path (where single subject EEG datasets and channel locations files are stored)
maindir = [pwd '\conn_decode_v1\' exp.settings '\']; %gets written over below
bdir = [maindir study_name '\'];

% Output directory (where decoding results will be saved)
% output_dir = [bdir 'SVR_Results\']; %for support vectore regression (SVR)
output_dir = [bdir 'SVR_Results_M\']; %for SVR with M (prob. of guess)

% if folder doesn't exist yet, create one
if ~exist(output_dir)
    mkdir(output_dir);
end

    
% Filepaths of single subject datasets (relative to the base directory)
sbj_code = cell(length(sbj_todo),1);
for ii = 1:length(sbj_todo)
    sbj_code{ii,1} = ['EEG_data\sbj' exp.participants{sbj_todo(ii)}];
end
clear ii

% Automatically calculates number of subjects from the number of data files
nsbj = size(sbj_code, 1);

% MATLAB workspace name for single subject data arrays and structures
data_struct_name = 'eeg_sorted_cond'; % Data arrays for use with DDTBOX must use this name as their MATLAB workspace variable name
  

%% EEG Dataset Information

nchannels = 465; % Number of channel combos
sampling_rate = 1000; % Data sampling rate in Hz

% Corresponds to the time of the event of interest (e.g. stimulus presentation) 
% relative to the start of the epoch (in ms) 
pointzero = 1000;

% For plotting single subject temporal decoding results 
% (not required if performing spatial or spatiotemporal decoding)
channel_names_file = 'channel_OrientExo_ML.mat'; % Name of the .mat file containing channel labels and channel locations
channellocs = maindir; % Path of the directory containing channel information file


%% Condition and Discrimination Group (dcg) Information

% Label each condition / category
% Usage: cond_labels{condition number} = 'Name of condition';
% Example: cond_labels{1} = 'Correct Responses';
% Example: cond_labels{2} = 'Error Responses';
% Condition label {X} corresponds to data in column X of the single subject
% data arrays.
% 1 = L_NI    2 = L_I    3 = R_NI    4 = R_I

cond_labels{1} = 'L_NI';
cond_labels{2} = 'L_I';
cond_labels{3} = 'R_NI';
cond_labels{4} = 'R_I';
       

% Discrimination groups
% Enter the condition numbers of the conditions used in classification analyses.
% If performing support vector regression, only one condition number is
% needed per dcg.
% SVR example: dcg{1} = [1]; to perform SVR on data from condition 1

dcg{1} = [1];
dcg{2} = [2]; 
dcg{3} = [3];
dcg{4} = [4]; 



% Support Vector Regression (SVR) condition labels
% Enter the array entry containing condition labels for each discrimination
% group number. The SVR_labels array contains multiple cells, each
% containing a list of SVR condition labels.
% Usage: svr_cond_labels{dcg} = [cell number in SVR_labels];
% Example: svr_cond_labels{1} = [2]; to use SVR labels in cell 2 for dcg 1

svr_cond_labels{1} = [1];
svr_cond_labels{2} = [2];
svr_cond_labels{3} = [3];
svr_cond_labels{4} = [4];
              


% Label each discrimination group
% Usage: dcg_labels{Discrimination group number} = 'Name of discrimination group'
% Example: dcg_labels{1} = 'Correct vs. Error Responses';
% 1 = L_NI    2 = L_I    3 = R_NI    4 = R_I

dcg_labels{1} = 'L_NI';
dcg_labels{2} = 'L_I';
dcg_labels{3} = 'R_NI';
dcg_labels{4} = 'R_I';



% This section automaticallly fills in various parameters related to dcgs and conditions 
ndcg = size(dcg, 2);
nclasses = size(dcg{1}, 2);      
ncond = size(cond_labels, 2);



%% Multivariate Classification/Regression Parameters

analysis_mode = 3; % ANALYSIS mode (1 = SVM classification with LIBSVM / 2 = SVM classification with LIBLINEAR / 3 = SVR with LIBSVM)
normalise_data = 1; % Normalise data for each feature prior to decoding? 1 = Yes / 0 = No
stmode = 3; % SPACETIME mode (1 = spatial / 2 = temporal / 3 = spatio-temporal)
avmode = 1; % AVERAGE mode (1 = no averaging; use single-trial data / 2 = use run-averaged data). Note: Single trials needed for SVR
zscore_convert = 0; % Convert data into z-scores before decoding? 0 = no / 1 = yes
cross_val_steps = 10; % How many cross-validation steps (if no runs available)?
n_rep_cross_val = 10; % How many repetitions of full cross-validation with re-ordered data?
perm_test = 1; % Run decoding using permuted condition labels? 0 = no / 1 = yes
permut_rep = 20; % How many repetitions of full cross-validation for permuted labels analysis?

%---------------------------------------------------------------------------
% Choose window size (10 = 12 steps)
% step size in feat calculation = 5 ms
window_width_ms = 10; % Width of sliding analysis window in ms
step_width_ms = 10; % Step size with which sliding analysis window is moved through the trial

% Feature weights extraction
feat_weights_mode = 1; % Extract feature weights? 0 = no / 1 = yes

% Single subject decoding results plotting
display_on = 1; % Display single subject decoding performance results? 0 = no / 1 = yes
perm_disp = 1; % Display the permuted labels decoding results in figure? 0 = no / 1 = yes
plotting_mode = 'classic'; % Plotting style. Current options are 'cooper' and 'classic'
x_tick_spacing_steps = step_width_ms; % Number of time steps between X axis time labels. If set to empty ([]) then plotting defaults are used.

% 'quiet mode' option to suppress text output to the command line
quiet_mode = 2; 
% 1 = Allow all text output to command line
% 2 = Show only important warnings and analysis related info (makes decoding run faster)
% 3 = No text output



%% Copy All Settings Into the cfg Structure
% No user input required in this section

cfg.bdir = bdir;
cfg.output_dir = output_dir;
cfg.sbj_code = sbj_code;
cfg.nsbj = nsbj;
cfg.data_struct_name = data_struct_name;
cfg.nchannels = nchannels;
cfg.channel_names_file = channel_names_file;
cfg.channellocs = channellocs;
cfg.sampling_rate = sampling_rate;
cfg.pointzero = pointzero;
cfg.cond_labels = cond_labels;
cfg.dcg = dcg;
cfg.dcg_labels = dcg_labels;
cfg.svr_cond_labels = svr_cond_labels;
cfg.ndcg = ndcg;
cfg.nclasses = nclasses;
cfg.ncond = ncond;
cfg.study_name = study_name;
cfg.cross = cross;
cfg.analysis_mode = analysis_mode;
cfg.stmode = stmode;
cfg.avmode = avmode;
cfg.window_width_ms = window_width_ms;
cfg.step_width_ms = step_width_ms;
cfg.zscore_convert = zscore_convert;
cfg.perm_test = perm_test;
cfg.cross_val_steps = cross_val_steps;
cfg.n_rep_cross_val = n_rep_cross_val;
cfg.permut_rep = permut_rep;
cfg.feat_weights_mode = feat_weights_mode;
cfg.display_on = display_on;
cfg.perm_disp = perm_disp;
cfg.normalise_data = normalise_data;
cfg.plotting_mode = plotting_mode;
cfg.x_tick_spacing_steps = x_tick_spacing_steps;
cfg.quiet_mode = quiet_mode;



%% Run SVR Analyses For Specified Subjects and dcgs

for dcg_set = 1:length(dcgs_for_analyses)
    
    clear dcg_todo;
    dcg_todo = dcgs_for_analyses{dcg_set};
        
    for sbj = 1:length(sbj_todo)

        % Save subject and dcg numbers into the configuration settings
        % structure
        cfg.sbj = exp.participants{sbj_todo(sbj)};
        cfg.dcg_todo = dcg_todo;
        
        % Set subject-specific filepaths for opening and saving files
        cfg.data_open_name = [bdir, (sbj_code{sbj}), '.mat'];
        cfg.data_save_name = [bdir, (sbj_code{sbj}), '_data.mat'];
        
        if contains(output_dir,'_Results_M')
            cfg.regress_label_name = [bdir, sbj_code{sbj}, '_regress_sorted_dataM.mat']; % Filepath for regression labels file
        else
            cfg.regress_label_name = [bdir, sbj_code{sbj}, '_regress_sorted_data.mat']; % Filepath for regression labels file
        end

        % Run the decoding analyses
        decoding_erp_mine(cfg);

    end % of for sbj
    clear sbj
    
end % of for dcg_set
clear dcg_set










