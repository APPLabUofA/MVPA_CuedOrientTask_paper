

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
% study_name = 'PACoz';
% study_name = 'PACplv';
% study_name = 'PACtort';
study_name = 'PACglm';

% Base directory path (where single subject EEG datasets and channel locations files are stored)
maindir = [pwd '\decode_v2\' exp.settings '\']; %gets written over below
bdir = [maindir study_name '\'];

% Output directory (where decoding results will be saved)
output_dir = [bdir 'SVR_Results\']; %for support vectore regression (SVR)

% Which group-level statistical analysis method was used??
% 1 = Global null and population prevalence tests based on the minimum statistic
% 2 = Global null testing using paired-samples t tests
group_level_analysis_method = 2; 

% -------------------------------------------------------------------------
%% Appropriate labels

if strcmpi(study_name,'lowFreq_amp')
    % lowFreq_amp data
    filename{1} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVR_LIBSVM_DCGR_I.mat'];
    filename{2} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVR_LIBSVM_DCGR_NI.mat'];
    filename{3} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVR_LIBSVM_DCGL_I.mat'];
    filename{4} = [study_name '_GROUPRES_NSBJ28_win20_steps20_av1_st3_SVR_LIBSVM_DCGL_NI.mat'];
else  
    % All other feature
    filename{1} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVR_LIBSVM_DCGR_I.mat'];
    filename{2} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVR_LIBSVM_DCGR_NI.mat'];
    filename{3} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVR_LIBSVM_DCGL_I.mat'];
    filename{4} = [study_name '_GROUPRES_NSBJ28_win15_steps15_av1_st3_SVR_LIBSVM_DCGL_NI.mat'];
end

plotname{1} = [study_name ': Informative Cue, Right Target'];
plotname{2} = [study_name ': Non-Informative, Right Target'];
plotname{3} = [study_name ': Informative Cue, Left Target'];
plotname{4} = [study_name ': Non-Informative, Left Target'];

figname{1} = [study_name '_SVR_R_I.fig'];
figname{2} = [study_name '_SVR_R_NI.fig'];
figname{3} = [study_name '_SVR_L_I.fig'];
figname{4} = [study_name '_SVR_L_NI.fig'];

% -------------------------------------------------------------------------
%% Location to save plots

saveFig = [bdir 'Figures\GROUPRES_SVR\']; % set save directory of data set

% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end


% -------------------------------------------------------------------------
%% Loop through files

for ii = 1:length(filename)

    % Full filepath of group results file
    group_results_file = [output_dir filename{ii}];

    % Load the data file to get ANALYSIS parameters
    load(group_results_file);
    
    clear group_results_file

    % /////////////////////////////////////////////////////////////////////////
    %% Settings For Group Decoding Performance Results (Classification Accuracy / SVR Performance)

    ANALYSIS.permdisp = 1; % Display results from permuted labels analyses in the figure as separate line? 0 = no / 1 = yes
    ANALYSIS.disp.sign = 1; % Mark statistically significant steps in results figure? 0 = no / 1 = yes
    ANALYSIS.plot_robust = 0; % Choose estimate of location to plot. 0 = arithmetic mean / 1 = trimmed mean / 2 = median
    % Note: You can plot the mean even if the data were originally plotted with
    % the median or trimmed mean when performing group-level analyses.
    % If you originally plotted the trimmed mean, then you can also plot the median in this script.

    ANALYSIS.disp.temporal_decoding_colormap = 'jet'; % Colormap for temporal decoding results scalp maps

    % Figure position on the screen
    PLOT.FigPos = [100, 100, 800, 400];

    % Background colour for the figure
    PLOT.background_colour = [1, 1, 1]; % RGB values. Default [1, 1, 1] (white)

    % Figure title settings
    PLOT.TitleFontSize = 18;
    PLOT.TitleFontWeight = 'Bold'; % 'Normal' (Regular) or 'b' / 'Bold'
    PLOT.TitleString = plotname{ii};


    %% Y Axis label properties
    PLOT.ylabel.FontSize = 16;
    PLOT.ylabel.FontWeight = 'Bold'; % 'Normal' (Regular) or 'b' / 'Bold'
    PLOT.ylabel.Text = 'Fisher-transformed correlation coeff';

%     PLOT.Y_min = -0.5; % Y axis lower bound (Fisher-Z-transformed correlation coefficient)
%     PLOT.Y_max = 0.5; % Y axis upper bound (Fisher-Z-transformed correlation coefficient)
%     PLOT.Ysteps = 0.2; % Interval between Y axis labels/tick marks
    
    PLOT.Y_min = -0.2; % Y axis lower bound (Fisher-Z-transformed correlation coefficient)
    PLOT.Y_max = 0.2; % Y axis upper bound (Fisher-Z-transformed correlation coefficient)
    PLOT.Ysteps = 0.1; % Interval between Y axis labels/tick marks

    % Determine Y axis ticks
%     PLOT.Ytick = [(PLOT.Y_min + .1) : PLOT.Ysteps : (PLOT.Y_max - .1)];
    PLOT.Ytick = [(PLOT.Y_min + .1) : PLOT.Ysteps : (PLOT.Y_max)];


    %% X Axis label properties
    PLOT.xlabel.FontSize = 16;
    PLOT.xlabel.FontWeight = 'Bold'; % 'Normal' (Regular) or 'b' / 'Bold'
    PLOT.xlabel.Text = 'Time-steps [ms]'; % X axis label text

    % Size of X tick labels
    PLOT.XY_tick_labels_fontsize = 14;

    % Number of time steps per X axis tick label
    PLOT.x_tick_spacing = 100;

    ANALYSIS.pointzero = 1000; %was 1001

    % Change x-axis scale
    if strcmpi(study_name,'lowFreq_amp')
        ANALYSIS.xaxis_scale(2,1:end) = 0 : ANALYSIS.step_width_ms : 999;
    else
        ANALYSIS.xaxis_scale(2,1:end) = 125 : ANALYSIS.step_width_ms : 874;
    end
    
    PLOT.X_min = 0; % X axis lower bound (first time point)
    PLOT.X_max = ANALYSIS.pointzero; % Maximum value of X axis value. 
    PLOT.Xsteps = ANALYSIS.step_width_ms; % Automated calculations (no input required)

    % Determine X axis ticks
    % PLOT.Xtick = [ANALYSIS.xaxis_scale(1,1) : PLOT.x_tick_spacing : ANALYSIS.xaxis_scale(1, end)];
    PLOT.XtickLabel = (0 : PLOT.x_tick_spacing : ANALYSIS.pointzero) - ANALYSIS.pointzero; 
    % PLOT.XtickLabel = ANALYSIS.xaxis_scale(2, 1 : PLOT.x_tick_spacing : end) - ANALYSIS.pointzero; 
    PLOT.Xtick = 0 : PLOT.x_tick_spacing : ANALYSIS.pointzero;

    %% Decoding results properties
    PLOT.Res.Line = '-ks'; % Line colour and style
    PLOT.Res.LineColour = 'black';

    % Error bar plotting settings
    PLOT.Res.Error = 'black'; % Line colour and style
    PLOT.Res.ErrorLineWidth = 0.5;
    PLOT.Res.ErrorLine = 'none'; % Disables lines between error bars across steps

    PLOT.Res.LineWidth = 1.5; % Default 2
    PLOT.Res.MarkerEdgeColor = 'k'; % Default 'k' (black)
    PLOT.Res.MarkerFaceColor = 'w'; % Default 'w' (white)
    PLOT.Res.MarkerSize = 5; % Default 5

    % Error bar plotting
    PLOT.Res.Error = 'k'; % Line colour and style. Default 'k' (black)
    PLOT.Res.ErrorLineWidth = 0.5; % Default 0.5
    PLOT.Res.ErrorLine = 'none'; % Entering 'none' disables lines between error bars across steps


    %% Permutation results properties
    PLOT.PermRes.Line = '-ks'; % Line colour and style
    PLOT.PermRes.LineColour = 'blue';

    % Error bar plotting settings
    PLOT.PermRes.Error = 'blue'; % Line colour and style
    PLOT.PermRes.ErrorLineWidth = 0.5;
    PLOT.PermRes.ErrorLine = 'none'; % Disables lines between error bars across steps

    % Error bar plotting settings
    PLOT.PermRes.LineWidth = 1.5; % Default 2
    PLOT.PermRes.MarkerEdgeColor = 'b'; % Default 'b' (black)
    PLOT.PermRes.MarkerFaceColor = 'w'; % Default 'w' (white)
    PLOT.PermRes.MarkerSize = 5; % Default 5

    % Error bars for chance/permuted labels results
    PLOT.PermRes.Error = 'b'; % Line colour and style. Default 'b' (black)
    PLOT.PermRes.ErrorLineWidth = 0.5; % Default 0.5
    PLOT.PermRes.ErrorLine = 'none'; % Entering 'none' disables lines between error bars across steps


    %% Decoding Performance Plot Annotations

    % Define properties of line showing event onset
    PLOT.PointZero.Color = [0.5, 0.5, 0.5]; % Colour of line denoting event onset. Default [0.5, 0.5, 0.5] (gray)
    PLOT.PointZero.LineWidth = 3; % Width of line denoting event onset. Default 3
    % PLOT.PointZero.Point = find(ANALYSIS.data(3,:) == 1);
    PLOT.PointZero.Point = ANALYSIS.pointzero;

    if isempty(PLOT.PointZero.Point) %if event onset is outside plotted time window
        PLOT.PointZero.Point = 0; %set to 0
    end

    % Define properties of statistical significance markers
    PLOT.Sign.LineWidth = 10; % Default 10

    % Alpha level of shaded statistical significance regions. Range 0-1 (values approaching 1 are more
    % opaque colours)
    PLOT.Sign.FaceAlpha = .3;

    % Positions of statistical significance markers
    if ANALYSIS.analysis_mode ~= 3 % If not using SVR

        PLOT.Sign.LinePos = [PLOT.Y_min + 0.5, PLOT.Y_max - 0.5];

    elseif ANALYSIS.analysis_mode == 3 % If using SVR

        PLOT.Sign.LinePos = [PLOT.Y_min, PLOT.Y_max];

    end % of if ANALYSIS.analysis_mode


    %% Set Plotting Colourmaps

    PLOT.Sign.LineColor = 'reddishpurple';

    % Set colour maps for plotting. Function provided by Dr Patrick Cooper (thanks Patrick!)
    plot_colour_map = dd_make_colour_maps( ...
        PLOT.Res.LineColour, ...
        PLOT.PermRes.LineColour, ...
        PLOT.Sign.LineColor);


    %% Plot the information time-course for each analysis
    for ana = 1:size(ANALYSIS.RES.mean_subj_acc, 1)

        temp_data(1,:) = ANALYSIS.RES.mean_subj_acc(ana,:);
        temp_se(1,:) = ANALYSIS.RES.se_subj_acc(ana,:);

        temp_perm_data(1,:) = ANALYSIS.RES.mean_subj_perm_acc(ana,:);
        temp_perm_se(1,:) = ANALYSIS.RES.se_subj_perm_acc(ana,:);

        % Open new figure
        figure('color', PLOT.background_colour, 'Position', PLOT.FigPos);

        % Mark statistically significant time windows        
        if ANALYSIS.disp.sign == 1

            for step = 1:size(temp_data, 2)

                % Add background shading if statistically significant...
                if ANALYSIS.RES.h(ana, step) == 1

                    % Note, in order to show the effect, we need to slightly broaden the
                    % patch by one sample (otherwise it is invisble - hence the -1 +1 here)
                    xaxis_loc = ANALYSIS.xaxis_scale(2,step);
                    if strcmpi(study_name,'lowFreq_amp')
                        x_data = [xaxis_loc - 10, xaxis_loc - 10, xaxis_loc + 10, xaxis_loc + 10];
                    else
                        x_data = [xaxis_loc - 7.5, xaxis_loc - 7.5, xaxis_loc + 7.5, xaxis_loc + 7.5];
                    end
                    y_data = [PLOT.Y_min, PLOT.Y_max, PLOT.Y_max, PLOT.Y_min];

                    % Settings for significance marker patch object (used
                    % for adding shaded background colours at significant time
                    % steps)
                    sig_markers = patch('xdata', x_data, 'ydata', y_data);
                    sig_markers.FaceAlpha = PLOT.Sign.FaceAlpha;
                    sig_markers.FaceColor = plot_colour_map(3,:);
                    sig_markers.EdgeColor = 'none';

                    hold on;

                end % of if ANALYSIS.RES.h
                
            end % of for step
            clear step
            
        end % of if ANALYSIS.disp.sign


        % Plot group-averaged decoding accuracy (or median, trimmed mean etc.)
        plot(ANALYSIS.xaxis_scale(2,:), temp_data, PLOT.Res.Line, ...
            'Color', plot_colour_map(1, :), ...
            'LineWidth', PLOT.Res.LineWidth, ...
            'MarkerEdgeColor', PLOT.Res.MarkerEdgeColor,...
            'MarkerFaceColor', PLOT.Res.MarkerFaceColor, ...
            'MarkerSize', PLOT.Res.MarkerSize);

        hold on;

        % Plot error bars
        errorbar(ANALYSIS.xaxis_scale(2,:),temp_data, temp_se, PLOT.Res.Error, ...
            'linestyle', PLOT.Res.ErrorLine, ...
            'linewidth', PLOT.Res.ErrorLineWidth);

        hold on;

        % Plot permutation decoding results
        plot(ANALYSIS.xaxis_scale(2,:),temp_perm_data, PLOT.PermRes.Line, ...
                    'Color', plot_colour_map(2, :), ...
                    'LineWidth', PLOT.PermRes.LineWidth, ...
                    'MarkerEdgeColor', PLOT.PermRes.MarkerEdgeColor,...
                    'MarkerFaceColor', PLOT.PermRes.MarkerFaceColor, ...
                    'MarkerSize', PLOT.PermRes.MarkerSize);

        hold on;

        % Plot error bars
        errorbar(ANALYSIS.xaxis_scale(2,:),temp_perm_data, temp_perm_se, PLOT.PermRes.Error, ...
            'linestyle', PLOT.PermRes.ErrorLine, ...
            'linewidth', PLOT.PermRes.ErrorLineWidth);


        %% Define Axis Labels, Point Zero, Title

        % Axis limits
        axis([PLOT.X_min, PLOT.X_max, PLOT.Y_min, PLOT.Y_max]);

        % Axis Labels
        xlabel([PLOT.xlabel.Text], 'FontSize', PLOT.xlabel.FontSize, 'FontWeight', PLOT.xlabel.FontWeight);
        ylabel([PLOT.ylabel.Text], 'FontSize', PLOT.ylabel.FontSize, 'FontWeight', PLOT.ylabel.FontWeight);

        % Title
        if size(ANALYSIS.DCG,1) == 1 % If did not perform cross-decoding

            title(PLOT.TitleString,'FontSize', PLOT.TitleFontSize,...
                'FontWeight', PLOT.TitleFontWeight,'Interpreter','none');

        elseif size(ANALYSIS.DCG,1) == 2 % If performed cross-decoding

            title([PLOT.TitleString, ANALYSIS.DCG{1}, 'to', ANALYSIS.DCG{2}, ' N=',  num2str(ANALYSIS.nsbj)],...
                'FontSize', PLOT.TitleFontSize, 'FontWeight', PLOT.TitleFontWeight);

        end % of if size

        % Mark point zero (data was time-locked to this event)
        line([PLOT.PointZero.Point, PLOT.PointZero.Point], [PLOT.Y_max, PLOT.Y_min], ...
            'Color', PLOT.PointZero.Color,...
            'LineWidth', PLOT.PointZero.LineWidth);



        %% Define X/Y Ticks and Axis Labels

        % X and Y ticks properties
        set(gca,'Ytick', PLOT.Ytick, ...
            'Xtick', PLOT.Xtick, ...
            'fontsize', PLOT.XY_tick_labels_fontsize, ...
            'fontname','Arial');

        % Convert X axis ticks to strings in a cell. This is to avoid the
        % weird X axis tick mislocalisation bug that sometimes occurs in
        % MATLAB
        clear XTickLabel_Cell;

        XTickLabel_Cell = cell(length(PLOT.XtickLabel),1); %pre-allocate
        for x_tick_number = 1:length(PLOT.XtickLabel)

            XTickLabel_Cell{x_tick_number} = PLOT.XtickLabel(x_tick_number);

        end % of for x_tick_number

        set(gca, 'XTickLabel', XTickLabel_Cell);

        % clear temp-data
        clear temp_* 

        % Remove top and right borders and associated X/Y ticks
        box off;

        % Save figure
        savefig([saveFig figname{ii}])

    end
    clear ana ANALYSIS FW_ANALYSIS group_results_file plot_colour_map...
        x_tick_number XTickLabel_Cell


end %loop through files
clear ii PLOT plotname filename saveFig fignname









