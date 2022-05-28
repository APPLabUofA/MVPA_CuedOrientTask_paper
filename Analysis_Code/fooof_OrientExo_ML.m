
clear variables


%% Load processing settings
load('byTargets_ML_v1_Settings.mat');
% load('byCues_ML_v1_Settings.mat');


% /////////////////////////////////////////////////////////////////////////
%% Load previously processed data segments
% /////////////////////////////////////////////////////////////////////////
% note that some settings for loading data is specified in the function itself
[ALLEEG,EEG] = LoadProcData_OrientExo_ML(exp);

eeglab redraw %reopen eeglab gui

% /////////////////////////////////////////////////////////////////////////
%% Save loaded data (if not done so)
% /////////////////////////////////////////////////////////////////////////

data_out_L = cell(1,length(exp.participants));    %pre-allocate
data_out_R = cell(1,length(exp.participants));    %pre-allocate
i_part = 0; %counter of subjects
for ii = 1:2:length(exp.participants)*2 %two sets per subj
    
    i_part = i_part + 1;%counter of subjects
    
    % Extract segments by target side
    data_out_L{i_part}(:,:,:) = squeeze(ALLEEG(ii).data(:,:,:));
    data_out_R{i_part}(:,:,:) = squeeze(ALLEEG(ii+1).data(:,:,:)); %2nd set
        
end
clear i_part ii

% Set-up and save relevant information
times_out = ALLEEG(1).times; % Get times variable to save
chanlocs = ALLEEG(1).chanlocs; % Get times variable to save
s_rate = ALLEEG(1).srate; % Get times variable to save
%save with version for large files
save([pwd '/data_out_all_' exp.settings '.mat'],'data_out_L','data_out_R',...
    'times_out','chanlocs','s_rate','-v7.3')

% /////////////////////////////////////////////////////////////////////////
%% Load data if saved previously
load([pwd '/data_out_all_' exp.settings '.mat']);

% load specific EEG dataset to make EEGLab happy
EEG = pop_loadset('004_LT_byTargets_ML_v1.set');
eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Prepare Data for fooof
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_data\' exp.settings];

% if folder doesn't exist yet, create one
if ~exist(fooof_dir)
    mkdir(fooof_dir);
end

% Use only relevant times
timewin = [-1000 0]; %after onset of last cue to target onset
timelim = find(times_out>=timewin(1) & times_out<=timewin(2));


% Separate data by trial type
for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [fooof_dir '\' exp.participants{i_part} '\'];
    % if folder doesn't exist yet, create one
    if ~exist(tmpfile_loc)
        mkdir(tmpfile_loc);
    end
    
    % extract participant data
    tmp_L = data_out_L{i_part};
    tmp_R = data_out_R{i_part};
    
    % extract by electrode
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %for doing only a selection of electrodes
        
        % By left target & cue type
        tmp_L_v(:,:) = squeeze(tmp_L(i_elect,timelim,(valid{i_part}(position{i_part}==0)==1)));
        tmp_L_iv(:,:) = squeeze(tmp_L(i_elect,timelim,(valid{i_part}(position{i_part}==0)==0)));
        
        % Calculate a power spectrum with Welch's method
        [psd, freqs] = pwelch(tmp_L_v, 500, [], [], s_rate);
        % Save the power spectra out to mat files
        save([tmpfile_loc '\spect_Lv_' exp.elec_names{ii}], 'freqs', 'psd');
        clear psd freqs
        
        % Calculate a power spectrum with Welch's method
        [psd, freqs] = pwelch(tmp_L_iv, 500, [], [], s_rate);
        % Save the power spectra out to mat files
        save([tmpfile_loc '\spect_Liv_' exp.elec_names{ii}], 'freqs', 'psd');
        clear psd freqs
        
        %.....................................................................
        
        % By right target & cue type
        tmp_R_v(:,:) = squeeze(tmp_R(i_elect,timelim,(valid{i_part}(position{i_part}==1)==1)));
        tmp_R_iv(:,:) = squeeze(tmp_R(i_elect,timelim,(valid{i_part}(position{i_part}==1)==0)));
        
        % Calculate a power spectrum with Welch's method
        [psd, freqs] = pwelch(tmp_R_v, 500, [], [], s_rate);
        % Save the power spectra out to mat files
        save([tmpfile_loc '\spect_Rv_' exp.elec_names{ii}], 'freqs', 'psd');
        clear psd freqs
        
        % Calculate a power spectrum with Welch's method
        [psd, freqs] = pwelch(tmp_R_iv, 500, [], [], s_rate);
        % Save the power spectra out to mat files
        save([tmpfile_loc '\spect_Riv_' exp.elec_names{ii}], 'freqs', 'psd');
        clear psd freqs
        
        clear tmp_L_v tmp_L_iv tmp_R_v tmp_R_iv
    end
    clear ii i_elect tmp_L tmp_R tmpfile_loc
end
clear i_part timewin timelim


% Save participant ids
nametmp = cell2mat(exp.participants);
save([fooof_dir '\participants'], 'nametmp');
clear nametmp

% Save electrode labels
chantmp = exp.elec_names;
save([fooof_dir '\chlabel'], 'chantmp');
clear chantmp


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load results from fooof and save as .mat
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_data\' exp.settings];

% Remember current folder
prime_dir = pwd;

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [results_dir '\' exp.participants{i_part} '\'];
    
    cd(tmpfile_loc) % Go to folder with saved results
    list = cellstr(ls('*.json')); % Get list of result file names
    
    % Loop through file list
    for ifile = 1:length(list)
        fname = list{ifile};
        dat = importdata(fname);
        
        fooof_results = []; %pre-allocate
        for ind = 1:length(dat)
            fooof_results = [fooof_results, jsondecode(dat{ind})];
        end
        clear dat ind

        save(replace(fname,'.json','.mat'),'fooof_results') % save results
        clear fname fooof_results
    end
    clear ifile tmpfile_loc list

end
clear i_part

cd(prime_dir) %return to main directory



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load fooof band .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_data\' exp.settings];

% Remember current folder
prime_dir = pwd;

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [results_dir exp.participants{i_part} '\'];
    
    list = cellstr(ls([tmpfile_loc 'bands_*'])); % Get list of result file names
    list = list(~contains(list,'json')); % Only the .mat files
    
    for ielect = 1:length(exp.elec_names)
        
        elist = list(contains(list,['_' exp.elec_names{ielect}])); % Select files by electrode

        % Loop through electrode specific file list
        for ifile = 1:length(elist)
            load([tmpfile_loc elist{ifile}]) %load data
            
            % Create string with condition name
            cond_name = char(extractBetween(elist{ifile},'s_',['_' exp.elec_names{ielect}]));

            % Create variables w/cond names 
            eval(strcat('beta1_trialpk_', cond_name, '{i_part}(:,ielect)=squeeze(beta1(:,1));'))
            eval(strcat('alpha_trialpk_', cond_name, '{i_part}(:,ielect)=squeeze(alphas(:,1));'))
            eval(strcat('theta_trialpk_', cond_name, '{i_part}(:,ielect)=squeeze(thetas(:,1));'))
            
            beta1_trialpk{ifile,i_part}(:,ielect) = squeeze(beta1(:,1));
            alpha_trialpk{ifile,i_part}(:,ielect) = squeeze(alphas(:,1));
            theta_trialpk{ifile,i_part}(:,ielect) = squeeze(thetas(:,1));

            clear alphas thetas beta1 cond_name
        end
        clear ifile elist
    end
    clear tmpfile_loc list ielect
end
clear i_part





% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load fooof fit results .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_data\' exp.settings];

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [results_dir exp.participants{i_part} '\'];
    
    list = cellstr(ls([tmpfile_loc 'f_results_*'])); % Get list of result file names
    list = list(~contains(list,'json')); % Only the .mat files
    
    % Loop through file list
    for ifile = 1:length(list)
        load([tmpfile_loc list{ifile}]) %load data
        
        cond_name = char(extractBetween(list{ifile},'s_','.mat'));
        
        % Get electrode id number
        chan_name = char(extractAfter(cond_name,'_'));
        ielect = find(matches(exp.elec_names,chan_name));
        
        % Get condition label
        con_lab = char(extractBefore(cond_name,'_'));
        
        eval(strcat('elect_track_', con_lab, '{ielect,i_part}=chan_name;')) %make sure things match
        
        for itrial = 1:length(fooof_results)
        
            eval(strcat('aperiodic_', con_lab, '{ielect,i_part}(1:2,itrial)=fooof_results(itrial).aperiodic_params_;'))
            eval(strcat('r_square_', con_lab, '{ielect,i_part}(:,itrial)=fooof_results(itrial).r_squared_;'))
            eval(strcat('fit_error_', con_lab, '{ielect,i_part}(:,itrial)=fooof_results(itrial).error_;'))
        
        end
        clear itrial con_lab ielect cond_name chan_name ielect fooof_results
    end
    clear ifile tmpfile_loc list
end
clear i_part results_dir fooof_dir


save(['fooof_fits_' exp.settings],'aperiodic*','r_square*','fit_error*',...
    'elect_track*')


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ------------------------------------------------------------------------- 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Load processed data

load('byTargets_ML_v1_Settings.mat');

load(['fooof_fits_' exp.settings]);

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Permutation test - aperiodic offset
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Aperiodic parameters: [Offset, Exponent]


test_data = NaN(length(exp.brainelecs),1,2,2,length(exp.participants)); %pre-allocate 
for ipart = 1:length(exp.participants)

    for ii = 1:length(exp.brainelecs)
%         ielect = R_Chan(ii);

        test_data(ii,1,2,1,ipart) = squeeze(median(aperiodic_Riv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,1,ipart) = squeeze(median(aperiodic_Liv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,2,2,ipart) = squeeze(median(aperiodic_Rv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,2,ipart) = squeeze(median(aperiodic_Lv{ii,ipart}(1,:),'omitnan'));

    end
    clear ielect ii

end
clear ipart 


%data  - 5D matrix of electrode x time points x side x cue x subjects
% dim 3 = [L R] 
% dim 4 = [no-info info]

alpha = 0.05/3; %Set alpha level
nperm = 1e4; %Number of permutations

% Aperiodic exponent
size(test_data)


dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


out = squeeze(mean(mean(test_data([7,11],1,:,:,:),3),5));

out = squeeze(std(mean(test_data([7,11],1,:,:,:),3),[],5));

clear test_data out






% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Permutation test - aperiodic exponent
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Aperiodic parameters: [Offset, Exponent]


test_data = NaN(length(exp.brainelecs),1,2,2,length(exp.participants)); %pre-allocate 
for ipart = 1:length(exp.participants)

    for ii = 1:length(exp.brainelecs)
%         ielect = R_Chan(ii);

        test_data(ii,1,2,1,ipart) = squeeze(median(aperiodic_Riv{ii,ipart}(2,:),'omitnan'));

        test_data(ii,1,1,1,ipart) = squeeze(median(aperiodic_Liv{ii,ipart}(2,:),'omitnan'));

        test_data(ii,1,2,2,ipart) = squeeze(median(aperiodic_Rv{ii,ipart}(2,:),'omitnan'));

        test_data(ii,1,1,2,ipart) = squeeze(median(aperiodic_Lv{ii,ipart}(2,:),'omitnan'));

    end
    clear ielect ii

end
clear ipart 


%data  - 5D matrix of electrode x time points x side x cue x subjects
% dim 3 = [L R] 
% dim 4 = [no-info info]

alpha = 0.05/3; %Set alpha level
nperm = 1e4; %Number of permutations

size(test_data)


dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


clear test_data out




% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Permutation test - r_square
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


test_data = NaN(length(exp.brainelecs),1,2,2,length(exp.participants)); %pre-allocate 
for ipart = 1:length(exp.participants)

    for ii = 1:length(exp.brainelecs)
%         ielect = R_Chan(ii);

        test_data(ii,1,2,1,ipart) = squeeze(median(r_square_Riv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,1,ipart) = squeeze(median(r_square_Liv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,2,2,ipart) = squeeze(median(r_square_Rv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,2,ipart) = squeeze(median(r_square_Lv{ii,ipart}(1,:),'omitnan'));

    end
    clear ielect ii

end
clear ipart 

% clear aperiodic* r_square* fit_error* elect_track*


%data  - 5D matrix of electrode x time points x side x cue x subjects
% dim 3 = [L R] 
% dim 4 = [no-info info]

alpha = 0.05/3; %Set alpha level
nperm = 1e4; %Number of permutations

size(test_data)


dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


clear test_data out



% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Permutation test - fit error
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////


test_data = NaN(length(exp.brainelecs),1,2,2,length(exp.participants)); %pre-allocate 
for ipart = 1:length(exp.participants)

    for ii = 1:length(exp.brainelecs)
%         ielect = R_Chan(ii);

        test_data(ii,1,2,1,ipart) = squeeze(median(fit_error_Riv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,1,ipart) = squeeze(median(fit_error_Liv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,2,2,ipart) = squeeze(median(fit_error_Rv{ii,ipart}(1,:),'omitnan'));

        test_data(ii,1,1,2,ipart) = squeeze(median(fit_error_Lv{ii,ipart}(1,:),'omitnan'));

    end
    clear ielect ii

end
clear ipart 

% clear aperiodic* r_square* fit_error* elect_track*


%data  - 5D matrix of electrode x time points x side x cue x subjects
% dim 3 = [L R] 
% dim 4 = [no-info info]

alpha = 0.05/4; %Set alpha level
nperm = 1e4; %Number of permutations

size(test_data)


dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);

clear dims test_results


clear test_data out nperm alpha


clear aperiodic* r_square* fit_error* elect_track*












