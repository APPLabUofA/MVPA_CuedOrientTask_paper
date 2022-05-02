function corrEEG = conn_corr_amp(eegData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%
%   >>  corrEEG = conn_corr(eegData);
%
% Input(s):
%           eegData  - input filtered dataset (channel x trials x time)
%   
% Outputs:
%           corrEEG  - output correlation matrix
%
% Info:
%           Computes Spearman's temporal correlation coefficient for each 
%           possible pair of EEG channels.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numChannels = size(eegData,1);
numTrials = size(eegData,2);

% get amplitude envelope
amp = NaN(numChannels,numTrials,size(eegData,3)); %pre-allocate
for kk = 1:numChannels
    amp(kk,:,:) = ft_preproc_hilbert(squeeze(eegData(kk,:,:)),'abs');
end
clear kk

corrEEG = NaN(numTrials,numChannels,numChannels); %pre-allocate
for ii = 1:numTrials
    corrEEG(ii,:,:) = corr(squeeze(amp(:,ii,:))','Type','Spearman');
end
clear ii




