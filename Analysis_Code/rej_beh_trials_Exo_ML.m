function [errordeg,valid,position] = rej_beh_trials_Exo_ML(exp,ALLEEG,BEH_path)
% Must load data using LoadProcData_OrientWheel_Exo_ML.m for function to work.
% ALLEEG structure should contain all rejected trials.
% Returns behavioral information on each trial with trials rejected during
% the EEG pre-processing pipeline removed. 

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Remove rejected trials
% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
errordeg = cell(length(exp.participants),1); %pre-allocate
valid = cell(length(exp.participants),1); %pre-allocate
position = cell(length(exp.participants),1); %pre-allocate

if contains(exp.setname,'byTargets','IgnoreCase',true)
    i_part = 0; %counter of subjects
else
    i_part = 1; %don't need counter for byCue settings
end

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Switch 
sets_name = split(exp.settings,'_'); %get what event trials are aligned to
switch sets_name{1}
% /////////////////////////////////////////////////////////////////////////
    
    case 'byTargets'
        
        for i_ALLEEG = 1:2:length(exp.participants)*2 %need every other file
            [n,m] = size(ALLEEG(i_ALLEEG).rejtrial);
            i_part = i_part + 1;
            %% Get list of rejected trials
            pip = 1;
            for ni = 1:n %for when there are more than 1 column
                for mi = 1:m
                    if ~isempty(ALLEEG(i_ALLEEG).rejtrial(ni,mi).ids)
                        rejlist{pip} = ALLEEG(i_ALLEEG).rejtrial(ni,mi).ids;
                        pip = 1 + pip;
                    end
                end
                clear mi
            end
            %% Load behavior data
            load([BEH_path num2str(exp.participants{i_part}) '_Orient_Exo.mat'])
%             ALLEEG(i_ALLEEG).error_deg

            %% Reject behavioral data trials
            if pip > 1 %if trials were rejected
                % Remove practice trials (first 20 trials)
                out_errordeg_temp = data.errorDegrees(1,21:end);
                out_valid_temp = data.valid(1,21:end);
                out_position_temp = data.position(1,21:end);
                
                % Need to remove trials from subject due to recording error
%                 if strcmpi(num2str(exp.participants{i_part}),'009')
%                     out_errordeg_temp(344:348) = [];
%                     out_valid_temp(344:348) = [];
%                     out_position_temp(344:348) = [];
%                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to remove trials from subject due to recording error
                if strcmpi(num2str(exp.participants{i_part}),'009')
                    out_errordeg_temp(344:end) = [];
                    out_valid_temp(344:end) = [];
                    out_position_temp(344:end) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Each set of rejected trials needs to be removed in order sequentially
                for mi = 1:length(rejlist)
                    tmplist = [rejlist{mi}];
                    out_errordeg_temp(tmplist) = []; %removes the trials
                    out_valid_temp(tmplist) = []; %removes the trials
                    out_position_temp(tmplist) = []; %removes the trials
                    clear tmplist
                end
                clear mi

            elseif pip == 1 %if no trials were rejected, rejlist variable not created
                % Remove practice trials (first 20 trials)
                out_errordeg_temp = data.errorDegrees(1,21:end);
                out_valid_temp = data.valid(1,21:end);
                out_position_temp = data.position(1,21:end);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to remove trials from subject due to recording error
                if strcmpi(num2str(exp.participants{i_part}),'009')
                    out_errordeg_temp(344:end) = [];
                    out_valid_temp(344:end) = [];
                    out_position_temp(344:end) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            % create variable with selected BEH 
            errordeg{i_part} = out_errordeg_temp;
            valid{i_part} = out_valid_temp;
            position{i_part} = out_position_temp;

            clear rejlist n m pip ni data prefs
            % clears variables that end/begin with...
            clear -regexp \<sgnrank_ _temp\>

        end
        clear i_part i_ALLEEG

% /////////////////////////////////////////////////////////////////////////    
%% by_Cues    
    otherwise 
        
        for i_ALLEEG = 1:length(exp.participants)
            [n,m] = size(ALLEEG(i_ALLEEG).rejtrial);
            i_part = i_ALLEEG;
            %% Get list of rejected trials
            pip = 1;
            for ni = 1:n %for when there are more than 1 column
                for mi = 1:m
                    if ~isempty(ALLEEG(i_ALLEEG).rejtrial(ni,mi).ids)
                        rejlist{pip} = ALLEEG(i_ALLEEG).rejtrial(ni,mi).ids;
                        pip = 1 + pip;
                    end
                end
                clear mi
            end
            %% Load behavior data
            load([BEH_path num2str(exp.participants{i_part}) '_Orient_Exo.mat'])

            %% Reject behavioral data trials
            if pip > 1 %if trials were rejected
                % Remove practice trials (first 20 trials)
                out_errordeg_temp = data.errorDegrees(1,21:end);
                out_valid_temp = data.valid(1,21:end);
                out_position_temp = data.position(1,21:end);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to remove trials from subject due to recording error
                if strcmpi(num2str(exp.participants{i_part}),'009')
%                     out_errordeg_temp(344:347) = [];
%                     out_valid_temp(344:347) = [];
%                     out_position_temp(344:347) = [];
                    out_errordeg_temp(344:end) = [];
                    out_valid_temp(344:end) = [];
                    out_position_temp(344:end) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Each set of rejected trials needs to be removed in order sequentially
                for mi = 1:length(rejlist)
                    tmplist = [rejlist{mi}];
                    out_errordeg_temp(tmplist) = []; %removes the trials
                    out_valid_temp(tmplist) = []; %removes the trials
                    out_position_temp(tmplist) = []; %removes the trials
                    clear tmplist
                end
                clear mi

            elseif pip == 1 %if no trials were rejected, rejlist variable not created
                % Remove practice trials (first 20 trials)
                out_errordeg_temp = data.errorDegrees(1,21:end);
                out_valid_temp = data.valid(1,21:end);
                out_position_temp = data.position(1,21:end);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to remove trials from subject due to recording error
                if strcmpi(num2str(exp.participants{i_part}),'009')
%                     out_errordeg_temp(344:348) = [];
%                     out_valid_temp(344:348) = [];
%                     out_position_temp(344:348) = [];
                    out_errordeg_temp(344:end) = [];
                    out_valid_temp(344:end) = [];
                    out_position_temp(344:end) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            % create variable with selected BEH 
            errordeg{i_part} = out_errordeg_temp;
            valid{i_part} = out_valid_temp;
            position{i_part} = out_position_temp;

            clear rejlist n m pip ni data prefs
            % clears variables that end/begin with...
            clear -regexp \<sgnrank_ _temp\>

        end
        clear i_part i_ALLEEG
% /////////////////////////////////////////////////////////////////////////
end %end switch
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
