% This is the group analysis configuration script for the DDTBOX group
% analysis tutorial.
% 
% This script calls decoding_erp.m and relies on the functions supplied in
% DDTBOX. The toolbox can be downloaded from https://github.com/DDTBOX/DDTBOX/releases


%% Housekeeping

% Clears the workspace and closes all figure windows
clear variables;
% close all;

%% Add necessary folders to path list
addpath(genpath('C:\Users\ssshe\Documents\MathLab\Personal_Folders\Sarah\Matlab_Functions\DDTBOX'))
% addpath(genpath('D:\ssshe\Documents\MathLab\Personal_Folders\Sarah\Matlab_Functions\DDTBOX'))

% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
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
sbj_list = [1:length(exp.participants)]; %not based on actual id #
% sbj_list = [1:2]; %not based on actual id #

% Filepaths of single subject datasets (relative to the base directory)
sbj_todo = strings(length(sbj_list),1);
for ii = 1:length(sbj_list)
    sbj_todo(ii,1) = string(exp.participants{sbj_list(ii)});
end
clear ii

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
% study_name = 'iCOH';
study_name = 'PLV';
% study_name = 'wPLI';
% study_name = 'PLI';

% Base directory path (where single subject EEG datasets and channel locations files are stored)
maindir = [pwd '\conn_decode_v1\' exp.settings '\']; %gets written over below
bdir = [maindir study_name '\'];

% Output directory (where decoding results will be saved)
output_dir = [bdir '\Decode_Results\']; %for classification 

% if folder doesn't exist yet, create one
if ~exist(output_dir)
    mkdir(output_dir);
end
   

%% EEG Dataset Information

nchannels = 465; % Number of channel combos
sampling_rate = 1000; % Data sampling rate in Hz

% Corresponds to the time of the event of interest (e.g. stimulus presentation) 
% relative to the start of the epoch (in ms) 
pointzero = 1000; 

channel_names_file = 'connchannel_OrientExo_ML.mat'; % Name of the .mat file containing channel labels and channel locations
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
       

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Discrimination groups
% Enter the condition numbers of the conditions used in classification analyses.
% Usage: dcg{discrimination group number} = [condition 1, condition 2];
% Example: dcg{1} = [1, 2]; to use conditions 1 and 2 for dcg 1

dcg{1} = [1, 2];
dcg{2} = [3, 4]; 

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Label each discrimination group
% Usage: dcg_labels{Discrimination group number} = 'Name of discrimination group'
% Example: dcg_labels{1} = 'Correct vs. Error Responses';
% 1 = L_NI    2 = L_I    3 = R_NI    4 = R_I

dcg_labels{1} = 'L_NI vs. L_I';
dcg_labels{2} = 'R_NI vs. R_I';


% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% This section automaticallly fills in various parameters related to dcgs and conditions 
ndcg = size(dcg, 2);
nclasses = size(dcg{1}, 2);      
ncond = size(cond_labels, 2);



%% Decoding Performance Analysis Parameters

% Specify the type of decoding analysis that was performed:
analysis_mode = 1; % ANALYSIS mode (1 = SVM classification with LIBSVM / 2 = SVM classification with LIBLINEAR / 3 = SVR with LIBSVM)
stmode = 3; % SPACETIME mode (1 = spatial / 2 = temporal / 3 = spatio-temporal)
avmode = 1; % AVERAGE mode (1 = no averaging; single-trial / 2 = run-averaged data) 

% Choose window size (10 = 12 steps)
% step size in feat calculation = 5 ms
window_width_ms = 10; % Width of sliding analysis window in ms
step_width_ms = 10; % Step size with which sliding analysis window is moved through the trial

% Specify alpha level
pstats = 0.05; % critical p-value

% Select group-level statistical analysis method 
% 1 = Global null and population prevalence tests based on the minimum statistic
% 2 = Global null testing using paired-samples t tests
% 3 = Global null testing (tmax) using paired-samples t tests for temporal data ONLY
group_level_analysis_method = 2; 

% In case set to 3 when not doing temporal analyses, change to 2
if stmode ~= 2 && group_level_analysis_method == 3
    group_level_analysis_method = 2;
end

% For spatial and spatiotemporal decoding results only:
% Last step (analysis time window) within the epoch to 
% include in group analyses. If left blank then this will
% be prompted at the command line while running group analyses.
% Example: laststep = [10] to perform tests on decoding results from analysis time windows 1 to 10.
laststep = 12; 

% For temporal decoding results only:
allchan = 1; % Are all possible channels analysed? 1 = yes (default for spatial or spatio-temporal decoding) / 2 = no
relchan = []; % Specify channels to be analysed (for temporal decoding only)


%__________________________________________________________________________
% If using the minimum statistic method for group-level analyses:

P2 = 1e7; % Number of second-level permutations to use; default: 100000 (recommend 1e7)

% Correct for multiple comparisons using the maximum statistic approach:
% 0 = no correction
% 1 = correction based on the maximum statistic (also applied to population prevalence estimates)
minstat_multcomp = 1; 


%__________________________________________________________________________
% If using paired-samples t tests for group-level analyses:

permstats = 2; % Testing against: 1 = theoretical chance level / 2 = permutation test results

% Testing against: 
% 1 = subject-averaged permuted labels decoding results (default)
% 2 = results of random permuted labels analysis repetitions drawn from each subject (stricter)
drawmode = 2; 

% Choose between two-tailed or one-tailed tests. 
% 'both' = two-tailed 
% 'right' or 'left' = one-tailed testing for above/below chance accuracy
groupstats_ttest_tail = 'right'; 

use_robust = 0; % Use Yuen's t, a robust version of the t test? 0 = no / 1 = yes
trimming = 20; % If using Yuen's t, select the trimming percentage for the trimmed mean

multcompstats = 3; % Correction for multiple comparisons: 
                    % 0 = no correction
                    % 1 = Bonferroni correction
                    % 2 = Holm-Bonferroni correction
                    % 3 = Strong FWER Control Permutation Test
                    % 4 = Cluster-Based Permutation Test
                    % 5 = KTMS Generalised FWER Control
                    % 6 = Benjamini-Hochberg FDR Control
                    % 7 = Benjamini-Krieger-Yekutieli FDR Control
                    % 8 = Benjamini-Yekutieli FDR Control
n_iterations = 1e4; %5000; % Number of permutation or bootstrap iterations for 
                      % resampling-based multiple comparisons correction procedures
ktms_u = 2; % u parameter for the KTMS GFWER control procedure
cluster_test_alpha = 0.05; % For cluster-based tests: Significance threshold 
                             % for inclusion of individual time windows into clusters




%% Decoding Performance Plotting Options

disp.on = 1; % Display a results figure? 0 = no / 1 = yes
permdisp = 1; % Display results from permuted labels analyses in the figure as separate line? 0 = no / 1 = yes
disp.sign = 1; % Mark statistically significant steps in results figure? 0 = no / 1 = yes
plot_robust = 0; % Choose estimate of location to plot. 0 = arithmetic mean / 1 = trimmed mean / 2 = median
plot_robust_trimming = 20; % Percentage trimming if using the trimmed mean
disp.temporal_decoding_colormap = 'jet'; % Colormap for temporal decoding scalp maps (default 'jet')

% Determine how many time steps between X axis ticks
% (e.g., with 12ms steps, a value of 5 means one X axis label every 60ms)
% (e.g., with 13ms steps, a value of 5 means one X axis label every 65ms)
x_tick_spacing_steps = 5; 

% Plotting style (options include 'cooper' and 'classic')
plotting_mode = 'classic';


%% Feature Weights Analysis Options

fw.do = 1; % Analyse feature weights? 0 = no / 1 = yes
fw.corrected = 1; % Use feature weights corrected using Haufe et al. (2014) method? 0 = no / 1 = yes

%---------------------------------------------------------------------------
% Time steps at which to perform statistical analyses on feature weights.
% Example: fw.steps_for_testing = [5:10]
% Input [] (empty vector) to manually input to the command line during FW analyses.
fw.steps_for_testing = [1:laststep]; 

    
%---------------------------------------------------------------------------
                           
fw.pstats = 0.05; % Alpha level for feature weights analyses
fw.use_robust = 0; % Use Yuen's t, a robust version of the t test? 0 = no / 1 = yes
fw.trimming = 20; % If using Yuen's t, select the trimming percentage for the trimmed mean
fw.ttest_tail = 'right'; % t test tail for feature weights analyses. 
                         %Should be set to 'right' for all standard analyses of FWs

fw.multcompstats = 3; % Feature weights correction for multiple comparisons:
                        % 1 = Bonferroni correction
                        % 2 = Holm-Bonferroni correction
                        % 3 = Strong FWER Control Permutation Test
                        % 4 = Cluster-Based Permutation Test (Currently not available)
                        % 5 = KTMS Generalised FWER Control
                        % 6 = Benjamini-Hochberg FDR Control
                        % 7 = Benjamini-Krieger-Yekutieli FDR Control
                        % 8 = Benjamini-Yekutieli FDR Control
fw.n_iterations = 1e4; %5000; % Number of permutation or bootstrap iterations for 
                        % resampling-based multiple comparisons correction procedures
fw.ktms_u = 0; % u parameter of the KTMS GFWER control procedure




%% Display Settings For Feature Weights Results 

fw.colormap = 'jet'; % Colormap for plotting of feature weights heat maps

% Consecutive time steps for which the feature weights matrix should be displayed
fw.disp_steps = [1:laststep];

% _________________________________________________________________________

% Display? 0 = no / 1 = yes
fw.display_matrix = 1; % Feature weights matrix

% _________________________________________________________________________
% Maps and stats averaged over selected analysis time windows
fw.display_average_zmap = 0; % Z-standardised absolute FWs
fw.display_average_uncorr_threshmap = 0; % Map of statistically significant FWs (uncorrected for multiple comparisons)
fw.display_average_corr_threshmap = 0; % Map of statistically significant FWs (corrected for multiple comparisons)

% _________________________________________________________________________
% Maps and stats for each selected analysis time window, plotted separately
fw.display_all_zmaps = 0; % Z-standardised absolute FWs
fw.display_all_uncorr_thresh_maps = 0; % Map of statistically significant FWs (uncorrected for multiple comparisons)
fw.display_all_corr_thresh_maps = 0; % Map of statistically significant FWs (corrected for multiple comparisons)

% _________________________________________________________________________
%% Copy All Settings Into ANALYSIS Structure

% This structure is passed as a single input argument to analyse_decoding_erp
% No user input is required for this section

ANALYSIS.bdir = bdir;
ANALYSIS.output_dir = output_dir;
ANALYSIS.nchannels = nchannels;
ANALYSIS.channel_names_file = channel_names_file;
ANALYSIS.channellocs = channellocs;
ANALYSIS.sampling_rate = sampling_rate;
ANALYSIS.pointzero = pointzero;
ANALYSIS.cond_labels = cond_labels;
ANALYSIS.dcg = dcg;
ANALYSIS.dcg_labels = dcg_labels;
ANALYSIS.ndcg = ndcg;
ANALYSIS.nclasses = nclasses;
ANALYSIS.ncond = ncond;
ANALYSIS.study_name = study_name;
ANALYSIS.sbjs_todo = sbj_todo;
% ANALYSIS.dcg_todo = dcg_todo;
ANALYSIS.cross = cross;
ANALYSIS.allchan = allchan;
ANALYSIS.relchan = relchan;
ANALYSIS.analysis_mode = analysis_mode;
ANALYSIS.stmode = stmode;
ANALYSIS.avmode = avmode;
ANALYSIS.window_width_ms = window_width_ms;
ANALYSIS.step_width_ms = step_width_ms;
ANALYSIS.laststep = laststep;
ANALYSIS.pstats = pstats;
ANALYSIS.group_level_analysis_method = group_level_analysis_method;
ANALYSIS.P2 = P2;
ANALYSIS.minstat_multcomp = minstat_multcomp;
ANALYSIS.permstats = permstats;
ANALYSIS.drawmode = drawmode;
ANALYSIS.groupstats_ttest_tail = groupstats_ttest_tail;
ANALYSIS.use_robust = use_robust;
ANALYSIS.trimming = trimming;
ANALYSIS.multcompstats = multcompstats;
ANALYSIS.n_iterations = n_iterations;
ANALYSIS.ktms_u = ktms_u;
ANALYSIS.cluster_test_alpha = cluster_test_alpha;
ANALYSIS.disp.on = disp.on;
ANALYSIS.permdisp = permdisp;
ANALYSIS.disp.sign = disp.sign;
ANALYSIS.disp.temporal_decoding_colormap = disp.temporal_decoding_colormap;
ANALYSIS.plot_robust = plot_robust;
ANALYSIS.plot_robust_trimming = plot_robust_trimming;
ANALYSIS.fw.do = fw.do;
ANALYSIS.fw.corrected = fw.corrected;
ANALYSIS.fw.steps_for_testing = fw.steps_for_testing;
ANALYSIS.fw.pstats = fw.pstats;
ANALYSIS.fw.use_robust = fw.use_robust;
ANALYSIS.fw.trimming = fw.trimming;
ANALYSIS.fw.ttest_tail = fw.ttest_tail;
ANALYSIS.fw.multcompstats = fw.multcompstats;
ANALYSIS.fw.n_iterations = fw.n_iterations;
ANALYSIS.fw.ktms_u = fw.ktms_u;
ANALYSIS.fw.display_matrix = fw.display_matrix;
ANALYSIS.fw.disp_steps = fw.disp_steps;
ANALYSIS.fw.colormap = fw.colormap;
ANALYSIS.fw.display_average_zmap = fw.display_average_zmap;
ANALYSIS.fw.display_average_uncorr_threshmap = fw.display_average_uncorr_threshmap;
ANALYSIS.fw.display_average_corr_threshmap = fw.display_average_corr_threshmap;
ANALYSIS.fw.display_all_zmaps = fw.display_all_zmaps;
ANALYSIS.fw.display_all_uncorr_thresh_maps = fw.display_all_uncorr_thresh_maps;
ANALYSIS.fw.display_all_corr_thresh_maps = fw.display_all_corr_thresh_maps;

ANALYSIS.disp.x_tick_spacing_steps = x_tick_spacing_steps;
ANALYSIS.disp.plotting_mode = plotting_mode;


%% Analyse Decoding Results For Specified Subjects and dcgs

for dcg_set = 1:length(dcg)
    
    if ANALYSIS.cross == 0
        ANALYSIS.dcg_todo = dcg_set;
    elseif ANALYSIS.cross == 1
        ANALYSIS.dcg_todo = dcgs_for_analyses{dcg_set};
    end

    analyse_decoding_erp_conn(ANALYSIS);

end
clear dcg_set










