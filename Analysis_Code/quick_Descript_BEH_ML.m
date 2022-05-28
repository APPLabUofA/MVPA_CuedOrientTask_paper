

%% Load processing settings
load('byTargets_ML_v1_Settings.mat');

%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);

%% Load raw behavioral data
for i_part = 1:length(exp.participants)
    
    % Load each subject's data
%     load(['M:\Data\OrientWheel_Exo\BEH\' exp.participants{i_part} '_Orient_Exo.mat'])
    load(['C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\BEH\' exp.participants{i_part} '_Orient_Exo.mat']) %my laptop
%     load(['C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\BEH\Staircasing\' exp.participants{i_part} '_stair_exo.mat']) %my laptop

    % Remove practice trials (first 20 trials)
    error_deg = data.errorDegrees(1,21:end);
    valid = data.valid(1,21:end);
    position = data.position(1,21:end);
%     targets = data.target(1,21:length(data.errorDegrees));
    
    errdeg_valid_L{i_part} = error_deg(valid == 1 & position == 0);
    errdeg_valid_R{i_part} = error_deg(valid == 1 & position == 1);
    errdeg_invalid_L{i_part} = error_deg(valid == 0 & position == 0);
    errdeg_invalid_R{i_part} = error_deg(valid == 0 & position == 1);
    
    % Record target color
    target_color(i_part) = prefs.targ_gray;

    clear error_deg valid data prefs position
end
clear i_part

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Find relative number of trials rejected during visual inspection

trl_count(1,:) = cell2mat({ALLEEG(1:end).trials}); %total trial count

vis_rej = cell(3,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
    vis_rej{1,ii} = ALLEEG(ii).rejtrial(3).ids;
    vis_rej{2,ii} = ALLEEG(ii).rejtrial(2).ids;
    vis_rej{3,ii} = ALLEEG(ii).rejtrial(1).ids;
end
clear ii

trl_count(2:4,:) = cellfun(@numel,vis_rej);

trl_final(1,:) = sum(trl_count(1:4,:));
trl_final(2,:) = sum(trl_count(2:4,:));
trl_final(3,:) = trl_count(2,:)./trl_final(1,:);
trl_final(4,:) = trl_count(2,:)./trl_final(2,:);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


for i_part = 1:length(exp.participants)
    
    % Load each subject's data
%     load(['M:\Data\OrientWheel_Exo\BEH\' exp.participants{i_part} '_Orient_Exo.mat'])
    load(['C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\BEH\' exp.participants{i_part} '_Orient_Exo.mat']) %my laptop

    % Remove practice trials (first 20 trials)
    error_deg = data.errorDegrees(1,21:end);
    valid = data.valid(1,21:end);
    position = data.position(1,21:end);
    
    errdeg_validcue_m(i_part) = mean(abs(error_deg(valid == 1)));
    errdeg_invalidcue_m(i_part) = mean(abs(error_deg(valid == 0)));
    
    errdeg_valid_L_m(i_part) = mean(abs(error_deg(valid == 1 & position == 0)));
    errdeg_valid_R_m(i_part) = mean(abs(error_deg(valid == 1 & position == 1)));
    errdeg_invalid_L_m(i_part) = mean(abs(error_deg(valid == 0 & position == 0)));
    errdeg_invalid_R_m(i_part) = mean(abs(error_deg(valid == 0 & position == 1)));
    
    errdeg_validcue{i_part} = error_deg(valid == 1);
    errdeg_invalidcue{i_part} = error_deg(valid == 0);
    
    errdeg_valid_L{i_part} = error_deg(valid == 1 & position == 0);
    errdeg_valid_R{i_part} = error_deg(valid == 1 & position == 1);
    errdeg_invalid_L{i_part} = error_deg(valid == 0 & position == 0);
    errdeg_invalid_R{i_part} = error_deg(valid == 0 & position == 1);

    clear error_deg valid data prefs position
end
clear i_part



mean(errdeg_validcue_m)
std(errdeg_validcue_m)

mean(errdeg_invalidcue_m)
std(errdeg_invalidcue_m)

[h,p]= ttest(errdeg_validcue_m,errdeg_invalidcue_m)

signrank(errdeg_validcue_m,errdeg_invalidcue_m)


cat_err_validcue = cat(2,errdeg_validcue{1:end});
cat_err_invalidcue = cat(2,errdeg_invalidcue{1:end});

figure
subplot(1,2,1); histogram(cat_err_validcue,20,'Normalization','probability')
ylim([0 0.35]); title('Valid Cue')
ylabel('Trial Proportion');xlabel('Degree Error')
subplot(1,2,2); histogram(cat_err_invalidcue,20,'Normalization','probability')
ylim([0 0.35]); title('Invalid Cue')
ylabel('Trial Proportion');xlabel('Degree Error')

figure
subplot(1,2,1); histogram(abs(cat_err_validcue),20,'Normalization','probability')
ylim([0 0.4]); 
ylabel('Trial Proportion');xlabel('Degree Error')
title('Valid Cue')
subplot(1,2,2); histogram(abs(cat_err_invalidcue),20,'Normalization','probability')
ylim([0 0.4]); 
title('Invalid Cue')
ylabel('Trial Proportion');xlabel('Degree Error')




mean(errdeg_valid_L_m)
std(errdeg_valid_L_m)

mean(errdeg_valid_R_m)
std(errdeg_valid_R_m)

mean(errdeg_invalid_L_m)
std(errdeg_invalid_L_m)

mean(errdeg_invalid_R_m)
std(errdeg_invalid_R_m)


[h,p]= ttest(errdeg_valid_L_m,errdeg_valid_R_m)
signrank(errdeg_valid_L_m,errdeg_valid_R_m)

[h,p]= ttest(errdeg_invalid_L_m,errdeg_invalid_R_m)
signrank(errdeg_invalid_L_m,errdeg_invalid_R_m)


[h,p]= ttest(errdeg_valid_L_m,errdeg_invalid_L_m)
signrank(errdeg_valid_L_m,errdeg_invalid_L_m)

[h,p]= ttest(errdeg_valid_R_m,errdeg_invalid_R_m)
signrank(errdeg_valid_R_m,errdeg_invalid_R_m)




cat_err_valid_L = cat(2,errdeg_valid_L{1:end});
cat_err_valid_R = cat(2,errdeg_valid_R{1:end});
cat_err_invalid_L = cat(2,errdeg_invalid_L{1:end});
cat_err_invalid_R = cat(2,errdeg_invalid_R{1:end});

figure
subplot(2,2,1); histogram(cat_err_valid_L,20,'Normalization','probability')
ylim([0 0.35]); 
title('Valid Cue: Left Target')
ylabel('Trial Proportion')
subplot(2,2,2); histogram(cat_err_valid_R,20,'Normalization','probability')
ylim([0 0.35]); 
title('Valid Cue: Right Target')
subplot(2,2,3); histogram(cat_err_invalid_L,20,'Normalization','probability')
ylim([0 0.35]); 
ylabel('Trial Proportion');xlabel('Degree Error')
title('Invalid Cue: Left Target')
subplot(2,2,4); histogram(cat_err_invalid_R,20,'Normalization','probability')
ylim([0 0.35]); 
title('Invalid Cue: Right Target')
xlabel('Degree Error')

figure
subplot(2,2,1); histogram(abs(cat_err_valid_L),20,'Normalization','probability')
ylim([0 0.5]); 
ylabel('Trial Proportion');
title('Valid Cue: Left Target')
subplot(2,2,2); histogram(abs(cat_err_valid_R),20,'Normalization','probability')
ylim([0 0.5]); 
title('Valid Cue: Right Target')
subplot(2,2,3); histogram(abs(cat_err_invalid_L),20,'Normalization','probability')
ylim([0 0.5]); 
ylabel('Trial Proportion');xlabel('Degree Error')
title('Invalid Cue: Left Target')
subplot(2,2,4); histogram(abs(cat_err_invalid_R),20,'Normalization','probability')
ylim([0 0.5]); 
xlabel('Degree Error')
title('Invalid Cue: Right Target')




cat_err_valid_L = cat(2,errdeg_valid_L{1:end});
cat_err_valid_R = cat(2,errdeg_valid_R{1:end});
cat_err_invalid_L = cat(2,errdeg_invalid_L{1:end});
cat_err_invalid_R = cat(2,errdeg_invalid_R{1:end});



