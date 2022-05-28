function FWanal_ConnDecoding_plot(ANALYSIS,FW,saveFig,figname,plotname,file_n,chanlocs)

%selecting the correct label
ii = file_n;


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
    channel_labels = cell(size(FW.AVERAGE_Z_DISP(:,:), 2),1); %pre-allocate
    for channel = 1:size(FW.AVERAGE_Z_DISP(:,:), 2)

        channel_labels{channel} = char(FW.chanlocs(1, channel).labels);

    end % of for channel

    % Channels plotted as rows, time-windows as columns
    resorted_data = [];
    resorted_data(:,:) = FW.AVERAGE_Z_DISP(:,:);
    resorted_data = resorted_data';

    % Create figure
    figure;
    colormap default
    clims = [-1 1]; %scaling of plot
    imagesc(resorted_data(:,:), clims);
%         imagesc(resorted_data(:,:));
    hold on;

    set(gca, 'Ytick', [1:size(FW.AVERAGE_Z_DISP, 2)]);
    set(gca, 'YTickLabel', (channel_labels));
    ylabel('Channel', 'FontSize', 5, 'FontWeight', 'b');

    set(gca, 'Xtick', [1:size(FW.AVERAGE_Z_DISP, 1)]);
    set(gca, 'XTickLabel', (FW.disp_steps));
    xlabel('Analysis time-step', 'FontSize', 9, 'FontWeight', 'b');

    % Color bar labels
    t = colorbar('peer',gca,'Ticks',[-1,-.5,0,.5,1]);
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
    ds.connectStrength = FW.AVERAGESTEPS_SELECT_FW_Z_MEAN';

    % replace very non-significant values with NaN
    ns = find(-0.45 < ds.connectStrength & ds.connectStrength < 0.45);
    ds.connectStrength(ns) = 0;

    ds.chanPairs = NaN(length(FW.chanlocs),2);
    for idx = 1:length(FW.chanlocs)
        ds.chanPairs(idx,:) = FW.chanlocs(idx).index;
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

    ds.connectStrength = FW.h_matrix_z_averagestep_corr;

    ds.chanPairs = NaN(length(FW.chanlocs),2);
    for idx = 1:length(FW.chanlocs)
        ds.chanPairs(idx,:) = FW.chanlocs(idx).index;
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

    for steps = 1:size(FW.p_matrix_z_corr, 2)

        % Replace non-brain electrodes with NaN
        to_plot = NaN(length(FW.chanlocs),1);
        to_plot(2:32) = FW.AVERAGE_Z_HEATS(steps, :)';

        figure;
        topoplot_decoding(to_plot, ...
            FW.chanlocs, 'style', 'both', 'electrodes', 'ptslabels', ...
            'maplimits', 'minmax', 'chaninfo', FW.chaninfo, ...
            'colormap', ANALYSIS.fw.colormap);

        hold on;

        % Color bar labels
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Z-standardised absolute feature weight');
        caxis([-1 1])
        t.Ticks = [-1,-.5,0,.5,1];

        title([plotname{ii} ': time step ', num2str(FW.steps_for_testing(steps))],...
            'FontSize', 8, 'FontWeight', 'b', 'Interpreter', 'none');

        saveas(gcf, [saveFig 'tstep_' num2str(FW.steps_for_testing(steps)) figname{ii}],...
            'fig')

        clear to_plot;

    end % of for steps

end % of if ANALYSIS.fw.display_all_zmaps
clear t

%% Display Heat Maps For Each Selected Step: Significance-Thresholded Maps Corrected For MCs

if ANALYSIS.fw.display_all_corr_thresh_maps == 1

    for steps = 1:size(FW.h_matrix_z_corr, 2)

        % Replace non-brain electrodes with NaN
        to_plot = NaN(length(FW.chanlocs),1);
        to_plot(2:32) = FW.h_matrix_z_corr{steps};

        figure;
        topoplot_decoding(to_plot, ... 
            FW.chanlocs, 'style', 'fill', 'electrodes', 'ptslabels', ...
            'numcontour', 1, 'conv', 'off', 'maplimits', [0 2], 'ccolor', [0 0 0], ...
            'ecolor', [1 1 1], 'chaninfo', FW.chaninfo, ...
            'colormap', ANALYSIS.fw.colormap_sig);

        hold on;

        % Color bar labels
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Feature weights corrected for MCs');

        title({['Feature weights MC corrected for time-step ',...
            num2str(FW.steps_for_testing(steps))]; plotname{ii}},...
            'FontSize', 8, 'FontWeight', 'b', 'Interpreter', 'none');

        % Save figure
        savefig([saveFig 'sigtopo_tstep_' num2str(FW.steps_for_testing(steps)) figname{ii}])

        clear to_plot;

    end % of for steps

end % of if ANALYSIS.fw.display_all_corr_thresh_maps
clear t












