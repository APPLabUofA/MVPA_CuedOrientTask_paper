function [plv] = conn_PLV(eegData)
% Computes the Phase Locking Value (PLV) for an EEG dataset.
%
% Input parameters:
%   eegData is a 3D matrix numChannels x numTrials x numTimePoints
%
% Output parameters:
%   plv is a 3D matrix - 
%     numTrials x numChannels x numChannels
%
%--------------------------------------------------------------------------
% This is to be done over a single time window instead of trials. See Chp
% 26 in Analyzing Neural Time Series Data by Cohen for more info.
% 
% NOTE:
% In order to extract the PLV between channels 17 and 20, use plv(:, 17, 20) 
% and NOT plv(:, 20, 17). The smaller channel number is to be used first.
%--------------------------------------------------------------------------
% 
% Reference:
%   Lachaux, J P, E Rodriguez, J Martinerie, and F J Varela. 
%   Measuring phase synchrony in brain signals. 
%   Human brain mapping 8, no. 4 (January 1999): 194-208. 
%   http://www.ncbi.nlm.nih.gov/pubmed/10619414.
% 
%--------------------------------------------------------------------------
% Written by: 
% Praneeth Namburi
% Cognitive Neuroscience Lab, DUKE-NUS
% 01 Dec 2009
% 
% Present address: Neuroscience Graduate Program, MIT
% email:           praneeth@mit.edu
% ADAPTED FROM THE ORIGINAL FUNCTION FOR FC LAB
% 
% Modified by SSS (Oct 2021)

numChannels = size(eegData, 1);
numPoint = size(eegData, 3);

% get phase
phaseData = NaN(numChannels,size(eegData,2),numPoint); %pre-allocate
for kk = 1:numChannels
    phaseData(kk,:,:) = ft_preproc_hilbert(squeeze(eegData(kk,:,:)),'angle');
end
clear kk

plv = NaN(size(eegData, 2), numChannels, numChannels); %pre-allocate
for channelCount = 1:numChannels-1
    channelData = squeeze(phaseData(channelCount, :, :));
    for compareChannelCount = channelCount+1:numChannels
        compareChannelData = squeeze(phaseData(compareChannelCount, :, :));
        
        plv(:, channelCount, compareChannelCount) = ...
            abs(sum(exp(1i*(channelData - compareChannelData)), 2))/numPoint;
        
        clear compareChannelData
    end
    clear compareChannelCount channelData
end
clear channelCount


plv = squeeze(plv);


