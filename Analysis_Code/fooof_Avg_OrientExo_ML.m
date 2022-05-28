

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
%% Prepare data for fooof
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_AvgData\' exp.settings];

% if folder doesn't exist yet, create one
if ~exist(fooof_dir)
    mkdir(fooof_dir);
end

% /////////////////////////////////////////////////////////////////////////

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
   
        
    % By left target & cue type
    tmp_L_v(:,:) = squeeze(mean(tmp_L(:,timelim,(valid{i_part}(position{i_part}==0)==1)),3));
    tmp_L_iv(:,:) = squeeze(mean(tmp_L(:,timelim,(valid{i_part}(position{i_part}==0)==0)),3));

    % Calculate a power spectrum with Welch's method
    [psd, freqs] = pwelch(tmp_L_v', 500, [], [], s_rate);
    % Save the power spectra out to mat files
    save([tmpfile_loc '\spect_Lv'], 'freqs', 'psd');
    clear psd freqs

    % Calculate a power spectrum with Welch's method
    [psd, freqs] = pwelch(tmp_L_iv', 500, [], [], s_rate);
    % Save the power spectra out to mat files
    save([tmpfile_loc '\spect_Liv'], 'freqs', 'psd');
    clear psd freqs

    %.....................................................................

    % By right target & cue type
    tmp_R_v(:,:) = squeeze(mean(tmp_R(:,timelim,(valid{i_part}(position{i_part}==1)==1)),3));
    tmp_R_iv(:,:) = squeeze(mean(tmp_R(:,timelim,(valid{i_part}(position{i_part}==1)==0)),3));

    % Calculate a power spectrum with Welch's method
    [psd, freqs] = pwelch(tmp_R_v', 500, [], [], s_rate);
    % Save the power spectra out to mat files
    save([tmpfile_loc '\spect_Rv'], 'freqs', 'psd');
    clear psd freqs

    % Calculate a power spectrum with Welch's method
    [psd, freqs] = pwelch(tmp_R_iv', 500, [], [], s_rate);
    % Save the power spectra out to mat files
    save([tmpfile_loc '\spect_Riv'], 'freqs', 'psd');
    clear psd freqs

    clear tmp_L_v tmp_L_iv tmp_R_v tmp_R_iv tmp_L tmp_R
end
clear i_part timewin timelim


% Save participant ids
nametmp = cell2mat(exp.participants);
save([fooof_dir '\participants'], 'nametmp');
clear nametmp




% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load results from fooof & save as .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_AvgData\' exp.settings];

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
%% Load fooof fit results .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_AvgData\' exp.settings];

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [results_dir '\' exp.participants{i_part} '\'];
    
    list = cellstr(ls([tmpfile_loc 'f_results_*'])); % Get list of result file names
    list = list(~contains(list,'json')); % Only the .mat files
    
    % Loop through file list
    for ifile = 1:length(list)
        load([tmpfile_loc list{ifile}]) %load data
        
        cond_name = char(extractBetween(list{ifile},'s_','.mat'));
        
        % Loop through brain electrodes
        aperiodic = NaN(length(exp.brainelecs),2); %pre-allocate
        r_square = NaN(length(exp.brainelecs),1); %pre-allocate
        fit_error = NaN(length(exp.brainelecs),1); %pre-allocate
        for ii = 1:length(exp.brainelecs)

            aperiodic(ii,1:2) = fooof_results(exp.brainelecs(ii)).aperiodic_params_;
            r_square(ii,1) = fooof_results(exp.brainelecs(ii)).r_squared_;
            fit_error(ii,1) = fooof_results(exp.brainelecs(ii)).error_;
            
        end
        clear ii fooof_results
        
        eval(strcat('aperiodic_', cond_name, '(:,:,i_part)=aperiodic;'))
        eval(strcat('r_square_', cond_name, '(:,:,i_part)=r_square;'))
        eval(strcat('fit_error_', cond_name, '(:,:,i_part)=fit_error;'))
        
        clear aperiodic r_square fit_error cond_name
    end
    clear ifile tmpfile_loc list
end
clear i_part


% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Permutation test
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

alpha = 0.05; %Set alpha level
nperm = 1e4; %Number of permutations

% Aperiodic exponent
test_data = NaN(size(aperiodic_Riv,1),size(aperiodic_Riv,3),1,2,2);
test_data(:,:,1,1,1) = aperiodic_Riv(:,2,:);
test_data(:,:,1,2,1) = aperiodic_Liv(:,2,:);
test_data(:,:,1,1,2) = aperiodic_Rv(:,2,:);
test_data(:,:,1,2,2) = aperiodic_Lv(:,2,:);
test_data = permute(test_data,[1 3 4 5 2]);
size(test_data)

%data  - 4D matrix of electrode x time points x conditions x subjects
dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


clear test_data


% -------------------------------------------------------------------------
% r-squared
test_data = NaN(size(r_square_Riv,1),size(r_square_Riv,3),1,2,2);
test_data(:,:,1,1,1) = r_square_Riv(:,1,:);
test_data(:,:,1,2,1) = r_square_Liv(:,1,:);
test_data(:,:,1,1,2) = r_square_Rv(:,1,:);
test_data(:,:,1,2,2) = r_square_Lv(:,1,:);
test_data = permute(test_data,[1 3 4 5 2]);
size(test_data)


%data  - 4D matrix of electrode x time points x conditions x subjects
dims = [3, 4]; % interaction
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


dims = 3; % side
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


dims = 4; % cue
test_results = calc_Fmax(test_data,[],dims,nperm,alpha);
[h, crit_p, adj_p] = fdr_bh(test_results.p,alpha,'pdep','yes');

clear dims adj_p crit_p h test_results


squeeze(mean(test_data,5))

squeeze(mean(mean(mean(test_data,5),4),3))



clear test_data




% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load fooof band .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
fooof_dir = [pwd '\fooof_AvgData\' exp.settings];

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

for i_part = 1:length(exp.participants)
    
    % Create new folder for each participant
    tmpfile_loc = [results_dir exp.participants{i_part} '\'];
    
    list = cellstr(ls([tmpfile_loc 'bands_*'])); % Get list of result file names
    list = list(~contains(list,'json')); % Only the .mat files
    
    % Loop through file list
    for ifile = 1:length(list)
        load([tmpfile_loc list{ifile}]) %load data
        
        cond_name = char(extractBetween(list{ifile},'s_','.mat'));
        
        % Create variables w/cond names 
        eval(strcat('betalow_pk_', cond_name, '(:,i_part)=beta1(exp.brainelecs,1);'))
        eval(strcat('alpha_pk_', cond_name, '(:,i_part)=alphas(exp.brainelecs,1);'))
        eval(strcat('theta_pk_', cond_name, '(:,i_part)=thetas(exp.brainelecs,1);'))
        
        clear alphas thetas beta1 cond_name
    end
    clear ifile tmpfile_loc list
end
clear i_part















