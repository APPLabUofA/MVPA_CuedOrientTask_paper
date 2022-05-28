function [ALLEEG,EEG] = LoadProcData_OrientExo_ML(exp)

%clear workspace
% ccc

% Load processing settings
% load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');

% For my laptop
FolderLoc = pwd;
% exp.pathname = 'C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\'; %data location
exp.pathname = 'D:\MathLab\Data\OrientWheel_Exo\'; %path on external hard drive
exp.electrode_locs = [FolderLoc '\EOG-electrode-locs-32_orientwheel.ced'];


% -------------------------------------------------------------------------
anal.tf = 'off'; % if loading TF data
anal.singletrials = 'off'; % if loading single trial data
anal.segments = 'on'; % if loading epochs
anal.tfelecs = exp.brainelecs; %electrodes
anal.singletrialselecs = exp.singletrialselecs; %single trial electrodes
% -------------------------------------------------------------------------
anal.rm_rejtrial_BEH = 'off'; % return response error w/rejected trials removed & model fit 
% BEH_path = 'M:\Data\OrientWheel_Exo\BEH\';
BEH_path = [exp.pathname 'BEH\']; %my computer
% -------------------------------------------------------------------------

nparts = length(exp.participants); %number of subjects
nsets = length(exp.setname); %will always be 1 in this study

% -------------------------------------------------------------------------
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% -------------------------------------------------------------------------



% #########################################################################
% #########################################################################
%% Load the data
% The main loop loops through sets, then participants, then events.
for i_set = 1:nsets
    exp.setname(i_set)
    tic
    
%     for i_part = nparts
    for i_part = 1:nparts
        sprintf(['Loading Participant ' num2str(exp.participants{i_part}) '...' ])
        
        % number of events of interest
        nevents = length(exp.events(i_set,:));
        
        part_name = exp.participants{i_part}; 
        
        for i_event = 1:nevents
            
            filename = [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
            
            % --------------------------------------------------------------------------------------------------------------------
            % Load the Time frequency data, if needed.
            if strcmpi('on',anal.tf) == 1 % only load these variables if we are loading time-frequency data
                
                %The variable ersp will be a 6D variable: (participants x sets x events x electrodes x frequencies x timepoints).
                ersp(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'ersp'));
                itc(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'itc'));

                if i_part == 1 && i_set == 1 && i_event == 1 %load time and freq data
                    times = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'times'));
                    freqs = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'freqs'));
                end
                
            end
            % --------------------------------------------------------------------------------------------------------------------
    
            % --------------------------------------------------------------------------------------------------------------------
            % Load the EEGLAB datasets, if needed.
            if strcmpi('on',anal.segments) == 1 || strcmp('on',anal.singletrials) == 1 % only load these variables if we are loading either ERP or single trial data
                try
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                catch
                    WaitSecs(.5)
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            elseif strcmpi('on',anal.tf) == 1 % if we are loading time-frequency data only, then we just need one of these.
                if i_part == 1 && i_set == 1 && i_event == 1
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            end
            % --------------------------------------------------------------------------------------------------------------------
            
            % --------------------------------------------------------------------------------------------------------------------
            % Load the Single Trial complex values, if needed
            if strcmpi('on',anal.singletrials) == 1 % only load these variables if we are loading single trial data
                
%                 ntrigs = length(exp.events{i_set});
                
                % Loads the time values and freqs
                if i_part == 1 && i_set == 1 && i_event == 1
                    times = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'times'));
                    freqs = struct2array(load([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' filename '.mat'],'freqs'));
                end
                
                %This block finds the latency of each event listed in exp.events, and the trials it appeared on.
%                 for nperevent = 1:ntrigs
%                     for i_trial = 1:EEG.trials
%                         if any(strcmp(num2str(exp.events{i_set,i_event}(nperevent)),EEG.epoch(i_trial).eventtype)) == 1
%                             targlatency{i_set,i_event,i_part} = [targlatency{i_set,i_event,i_part} EEG.epoch(i_trial).eventlatency(find(strcmp(num2str(exp.events{i_set,i_event}(nperevent)),EEG.epoch(i_trial).eventtype)))];
%                             event_trials{i_set,i_event,i_part} = [event_trials{i_set,i_event,i_part} i_trial];
%                         end
%                     end
%                 end

                % Load single trial data for each electrode
                for ii = 1:length(exp.singletrialselecs)
                    i_chan = exp.singletrialselecs(ii);
                    
                    % all_ersp is (participant x electrode).trials(freq x time x trial)
                    try %Unfortunately, this load procedure can break sometimes in a non-reproducible way. So if an error happens here, we wait half a second and try again.
                        channeldata = load([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
                        all_ersp(i_event,i_part,i_chan) = struct2cell(channeldata);
                    catch
                        WaitSecs(.5)
                        channeldata = load([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
                        all_ersp(i_event,i_part,i_chan) = struct2cell(channeldata);
                    end
                    clear channeldata
                    
                 end
            % --------------------------------------------------------------------------------------------------------------------
            end
        end
    end
    toc
end
clear i_chan i_event i_set i_part filename nevents nparts nsets ii part_name


eeglab redraw


% #########################################################################
% #########################################################################
%% Return BEH data corrected for rejected trials
% Get degree response errors on trials with the trials rejected during the
% EEG pre-processing removed & model fit parameters
if strcmpi('on',anal.rm_rejtrial_BEH) == 1
    [errordeg,valid,position] = rej_beh_trials_Exo_ML(exp,ALLEEG,BEH_path);
    
    %save processed behavior data and epoch information
    chanlocs = EEG.chanlocs; %to save electrode location info
    save([FolderLoc '\BEH_' exp.settings '.mat'],'errordeg','valid',...
        'position','chanlocs')
    
end
clear BEH_path chanlocs

% #########################################################################
% #########################################################################


end



