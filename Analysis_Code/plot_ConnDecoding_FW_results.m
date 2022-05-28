
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
% study_name = 'corr';
% study_name = 'dwPLI';
% study_name = 'iCOH';
study_name = 'PLV';

stmode = 3; % SPACETIME mode (1 = spatial / 2 = temporal / 3 = spatio-temporal)

% Base directory path (where single subject EEG datasets and channel locations files are stored)
maindir = [pwd '\conn_decode_v1\' exp.settings '\']; %gets written over below
bdir = [maindir study_name '\'];

% Output directory (where decoding results will be saved)
output_dir = [bdir 'Decode_Results\']; %for classification 

% Load actual channel info for topoplots
load([maindir 'channel_OrientExo_ML.mat'],'chanlocs')
chanlocs(33:34) = []; %remove EOGs
chanlocs(1) = []; %remove M2

% Which group-level statistical analysis method was used??
% 1 = Global null and population prevalence tests based on the minimum statistic
% 2 = Global null testing using paired-samples t tests
group_level_analysis_method = 2; 

% -------------------------------------------------------------------------
%% Appropriate labels

filename{1} = [study_name '_GROUPRES_NSBJ28_win10_steps10_av1_st3_SVM_LIBSVM_DCGR_NI vs. R_I.mat'];
filename{2} = [study_name '_GROUPRES_NSBJ28_win10_steps10_av1_st3_SVM_LIBSVM_DCGL_NI vs. L_I.mat'];


% Names to save plots
figname{1} = [study_name '_st3_SVM_R_NIvsI.fig'];
figname{2} = [study_name '_st3_SVM_L_NIvsI.fig'];
    
plotname{1} = [study_name ': Non-Informative vs. Informative Cue, Right Target'];
plotname{2} = [study_name ': Non-Informative vs. Informative Cue, Left Target'];

% -------------------------------------------------------------------------
%% Location to save plots

saveFig = [bdir 'Figures\GROUPRES_SVM\']; % set save directory of data set

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
    
    % ---------------------------------------------------------------------
    %% Change time steps to actual ms
    ANALYSIS.xaxis_scale(2,1:end) = 187.5 : (ANALYSIS.step_width_ms*5) : 787;
    
    % --------------------------------------------------------------------
    %% Feature Weights Results Plotting Settings

    % Display? 0 = no / 1 = yes
    ANALYSIS.fw.display_matrix = 1; % Feature weights matrix
    
    % Maps and stats for averaged analysis time windows
    ANALYSIS.fw.display_average_zmap = 1; % z-standardised average FWs
    ANALYSIS.fw.display_average_corr_threshmap = 1; % thresholded map t-test results corrected for multiple comparisons

    % Maps and stats for each analysis time window
    ANALYSIS.fw.display_all_zmaps = 0; % z-standardised average FWs
    ANALYSIS.fw.display_all_corr_thresh_maps = 0; % thresholded map t-test results corrected for multiple comparisons

    % Extra plotting options:
%     ANALYSIS.fw.colormap = redblue(256); % Colormap for plotting of feature weights scalp maps
    ANALYSIS.fw.colormap = redblue(64); % Colormap for plotting of feature weights scalp maps
    ANALYSIS.fw.colormap_sig = 'jet'; % Colormap for plotting significance-thresholded maps corrected for MCs
    
    
    %% Display Matrix of All Steps: Z-Standardised Absolute Feature Weights
    % This matrix is plotted from
    % FW_ANALYSIS.AVERAGE_Z_DISP{analysis-time-steps, channel}. Note that this
    % can be different from the statistically tested analysis time-windows 

    if ANALYSIS.fw.display_matrix == 1
        
        % Create labels
        channel_labels = cell(size(FW_ANALYSIS.AVERAGE_Z_DISP(:,:), 2),1); %pre-allocate
        for channel = 1:size(FW_ANALYSIS.AVERAGE_Z_DISP(:,:), 2)

            channel_labels{channel} = char(FW_ANALYSIS.chanlocs(1, channel).labels);

        end % of for channel

        % Channels plotted as rows, time-windows as columns
        resorted_data = [];
        resorted_data(:,:) = FW_ANALYSIS.AVERAGE_Z_DISP(:,:);
        resorted_data = resorted_data';

        % Create figure
        figure;
        colormap default
        clims = [-1 1]; %scaling of plot
        imagesc(resorted_data(:,:), clims);
%         imagesc(resorted_data(:,:));
        hold on;

        set(gca, 'Ytick', [1:size(FW_ANALYSIS.AVERAGE_Z_DISP, 2)]);
        set(gca, 'YTickLabel', (channel_labels));
        ylabel('Channel', 'FontSize', 6, 'FontWeight', 'b');

        set(gca, 'Xtick', [1:size(FW_ANALYSIS.AVERAGE_Z_DISP, 1)]);
        set(gca, 'XTickLabel', (FW_ANALYSIS.disp_steps));
        xlabel('Analysis time-step', 'FontSize', 9, 'FontWeight', 'b');
        
        % Color bar labels
        t = colorbar('peer',gca,'Ticks',[-1,-.5,0,.5,1]);
%         t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Z-standardised absolute feature weights');

        title({'Z-standardised absolute feature weights'; plotname{ii}},...
            'FontSize', 10, 'FontWeight', 'b', 'Interpreter', 'none');
        
        % Save figure
        savefig([saveFig 'matrix_' figname{ii} ])

    end % of if ANALYSIS.fw.display_matrix
    clear clims t resorted_data
    
    %% Display Average Heat Map For Selected Steps: Z-Standardised Absolute Feature Weights

    if ANALYSIS.fw.display_average_zmap == 1

        ds.connectStrengthLimits = [-1,1]; %min and max values of connection strengths to display
        ds.connectStrength = FW_ANALYSIS.AVERAGESTEPS_SELECT_FW_Z_MEAN';
        
        % replace very non-significant values with NaN
        ns = find(-0.3 < ds.connectStrength & ds.connectStrength < 0.3);
        ds.connectStrength(ns) = 0;
        
        ds.chanPairs = NaN(length(FW_ANALYSIS.chanlocs),2);
        for idx = 1:length(FW_ANALYSIS.chanlocs)
            ds.chanPairs(idx,:) = FW_ANALYSIS.chanlocs(idx).index;
        end
        clear idx

        figure;
        topoplot_connect(ds,chanlocs,ANALYSIS.fw.colormap);
        
        hold on;
        
        % Color bar labels
        colormap(gcf,ANALYSIS.fw.colormap)
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Z-standardised absolute feature weights');
        t.Ticks = [-1,-.5,0,.5,1];
        caxis([-1 1])
            
        title(plotname{ii}, 'FontSize', 8, 'FontWeight', 'b', 'Interpreter', 'none');
        
        saveas(gcf, [saveFig 'Ztopoplot_' figname{ii}], 'fig')

        clear ds

    end % of if ANALYSIS.fw.display_average_zmap
    clear t
    
    %% Display Average Heat Map For Selected Steps: Significance-Thresholded Map Corrected For MCs

    if ANALYSIS.fw.display_average_corr_threshmap == 1

        ds.connectStrength = FW_ANALYSIS.h_matrix_z_averagestep_corr;
        
        ds.chanPairs = NaN(length(FW_ANALYSIS.chanlocs),2);
        for idx = 1:length(FW_ANALYSIS.chanlocs)
            ds.chanPairs(idx,:) = FW_ANALYSIS.chanlocs(idx).index;
        end
        clear idx

        figure;
        topoplot_connect(ds,chanlocs,ANALYSIS.fw.colormap_sig);

        hold on;
        
        % Color bar labels
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Feature weights corrected for MCs');
        
        title({'Feature weights corrected for MCs (averaged across time-steps)';...
            plotname{ii}}, 'FontSize',8,'FontWeight','b','Interpreter', 'none');
        
        % Save figure
        savefig([saveFig 'sigtopo_' figname{ii}])

        clear ds

    end % of if ANALYSIS.fw.display_average_corr_threshmap
    clear t

    %% Display Heat Map For Each Selected Step: Z-Standardised Absolute Feature Weights

    if ANALYSIS.fw.display_all_zmaps == 1

        for steps = 1:size(FW_ANALYSIS.p_matrix_z_corr, 2)

            % Replace non-brain electrodes with NaN
            to_plot = NaN(length(FW_ANALYSIS.chanlocs),1);
            to_plot(2:32) = FW_ANALYSIS.AVERAGE_Z_HEATS(steps, :)';
    %         to_plot = FW_ANALYSIS.AVERAGE_Z_HEATS(steps, :);
    %         to_plot = to_plot';

            figure;
            topoplot_decoding(to_plot, ...
                FW_ANALYSIS.chanlocs, 'style', 'both', 'electrodes', 'ptslabels', ...
                'maplimits', 'minmax', 'chaninfo', FW_ANALYSIS.chaninfo, ...
                'colormap', ANALYSIS.fw.colormap);

            hold on;
            
            % Color bar labels
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Z-standardised absolute feature weight');
            caxis([-1 1])
            t.Ticks = [-1,-.5,0,.5,1];
        
            title([plotname{ii} ': time step ', num2str(FW_ANALYSIS.steps_for_testing(steps))],...
                'FontSize', 8, 'FontWeight', 'b', 'Interpreter', 'none');
            
            saveas(gcf, [saveFig figname{ii} '_tstep_' num2str(FW_ANALYSIS.steps_for_testing(steps))],...
                'fig')

            clear to_plot;

        end % of for steps

    end % of if ANALYSIS.fw.display_all_zmaps
    clear t
    
    %% Display Heat Maps For Each Selected Step: Significance-Thresholded Maps Corrected For MCs

    if ANALYSIS.fw.display_all_corr_thresh_maps == 1

        for steps = 1:size(FW_ANALYSIS.h_matrix_z_corr, 2)

            % Replace non-brain electrodes with NaN
            to_plot = NaN(length(FW_ANALYSIS.chanlocs),1);
            to_plot(2:32) = FW_ANALYSIS.h_matrix_z_corr{steps};
    %         to_plot(:,:) = FW_ANALYSIS.h_matrix_z_corr{steps};

            figure;
            topoplot_decoding(to_plot, ... 
                FW_ANALYSIS.chanlocs, 'style', 'fill', 'electrodes', 'ptslabels', ...
                'numcontour', 1, 'conv', 'off', 'maplimits', [0 2], 'ccolor', [0 0 0], ...
                'ecolor', [1 1 1], 'chaninfo', FW_ANALYSIS.chaninfo, ...
                'colormap', ANALYSIS.fw.colormap_sig);

            hold on;
            
            % Color bar labels
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Feature weights corrected for MCs');
            
            title({['Feature weights MC corrected for time-step ',...
                num2str(FW_ANALYSIS.steps_for_testing(steps))]; plotname{ii}},...
                'FontSize', 8, 'FontWeight', 'b', 'Interpreter', 'none');
            
            % Save figure
            savefig([saveFig figname{ii} '_sigtopo_tstep_'...
                num2str(FW_ANALYSIS.steps_for_testing(steps))])

            clear to_plot;

        end % of for steps

    end % of if ANALYSIS.fw.display_all_corr_thresh_maps
    clear t
    
    
end % of file loop
clear ii plotname filename saveFig figname study_name


clear ANALYSIS FW_ANALYSIS
























