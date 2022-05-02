function [ANALYSIS] = tmax_classifier_accuracies(ANALYSIS)
% Runs mxt_perm1 --> One sample permutation test based on a t-statistic and 
% null hypothesis of a mean of zero.
% Comes from the Mass_Univariate_ERP_Toolbox:
% 
% Groppe, D.M., Urbach, T.P., Kutas, M. (2011) Mass univariate analysis of 
% event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology. 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x
% 
% 
% ANALYSIS.drawmode == 2 % against one randomly drawn value (from all cross-val 
% repetitions for each participant) for stricter test
% from drawing values from random-labels test average permutation results 
%    across cross-validation steps, but draw later one for each participant 
%    for statistical testing
% ANALYSIS.RES.all_subj_perm_acc_reps_draw{subj, ana/chan, step}

% ANALYSIS.RES.draw_subj_perm_acc(sbj, na, step)
% used in actual statistical testing


% Needs to be transformed
ANALYSIS.RES.mean_subj_acc = ANALYSIS.RES.mean_subj_acc';
ANALYSIS.RES.se_subj_acc = ANALYSIS.RES.se_subj_acc';
ANALYSIS.RES.mean_subj_perm_acc = ANALYSIS.RES.mean_subj_perm_acc';
ANALYSIS.RES.se_subj_perm_acc = ANALYSIS.RES.se_subj_perm_acc';

for na = 1:size(ANALYSIS.RES.mean_subj_acc, 1) % analysis/channel
    for step = 1:size(ANALYSIS.RES.mean_subj_acc, 2) % step/analysis time window
        for sbj = 1:ANALYSIS.nsbj

            temp = randperm(size(ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj, na, step}(:,:), 2));
            drawone = temp(1); clear temp;
            ANALYSIS.RES.draw_subj_perm_acc(sbj, na, step) = ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj, na, step}(1, drawone);
            clear drawone;

        end % of for sbj
        clear sbj
    end % of for step
    clear step
end % of for na/chan
clear na


%  data   - 3D matrix of data (Channel x Time x Participant)
real_decoding_scores = permute(ANALYSIS.RES.all_subj_acc(:, :, :),[2 3 1]); %rearrange dims
perm_decoding_scores = permute(ANALYSIS.RES.draw_subj_perm_acc(:, :, :),[2 3 1]); %rearrange dims

data = real_decoding_scores - perm_decoding_scores; %test difference vs 0

% Get type of ttest
if strcmpi('right',ANALYSIS.groupstats_ttest_tail)
    tail = 1; %upper tailed test
elseif strcmpi('left',ANALYSIS.groupstats_ttest_tail)
    tail = -1; %lower tailed test
elseif strcmpi('both',ANALYSIS.groupstats_ttest_tail)
    tail = 0; %two tailed test
end
    

% Run permutation tmax test
[pval, ~, tmx_ptile, corrected_p, corrected_h, ~, ~] = mxt_perm1_Phipson_Smyth(data, ANALYSIS.n_iterations, ANALYSIS.pstats, tail);

h_ttest = pval < ANALYSIS.pstats; %creates a logical for significant values


%% To be returned by function
ANALYSIS.RES.h_ttest(:, :) = h_ttest;
ANALYSIS.RES.p_ttest(:,:) = pval;
ANALYSIS.RES.critical_t = tmx_ptile;
ANALYSIS.RES.corrected_p(:, :) = corrected_p;
ANALYSIS.RES.corrected_h(:,:) = corrected_h;

% Marking h values (statistical significance) for plotting
ANALYSIS.RES.h = ANALYSIS.RES.h_ttest;







