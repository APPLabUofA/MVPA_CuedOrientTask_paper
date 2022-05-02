%% Imaginary Coherence
% 
% [icoh] = conn_iCOH(eegData) 
% 
% Input(s):
%           eegData  - input filtered dataset (channel x trials x time)
% 
% Outputs:
%           icoh     - imaginary coherence matrix (trials x channel x channel)
% 
%--------------------------------------------------------------------------
% NOTE:
% In order to extract the ppc between channels 17 and 20, use ppc(:,17,20) 
% and NOT ppc(:,20,17). The smaller channel number is to be used first.
%--------------------------------------------------------------------------

function [icoh] = conn_iCOH(eegData) 

numChannels = size(eegData,1);
numTrials = size(eegData,2);
numPoint = size(eegData,3);

% time needs to be first
sig = permute(eegData, [3 2 1]);

% get complex number
complex_sig = NaN(numChannels,numPoint,numTrials); %pre-allocate
for kk = 1:numChannels
    complex_sig(kk,:,:) = hilbert(squeeze(sig(:,:,kk)));
end
clear kk sig

icoh = NaN(numTrials, numChannels, numChannels); %pre-allocate
for channelCount = 1:numChannels-1
    channelData = squeeze(complex_sig(channelCount, :, :));
    for compareChannelCount = channelCount+1:numChannels
        compareChannelData = squeeze(complex_sig(compareChannelCount, :, :));
        
        % cross-spectral density
        spec1 = sum(channelData.*conj(channelData),1);
        spec2 = sum(compareChannelData.*conj(compareChannelData),1);
        specX = sum(channelData.*conj(compareChannelData),1);
        
        % imaginary coherence
        icoh(:, channelCount, compareChannelCount) = abs(imag(specX./sqrt(spec1.*spec2)));
        
        clear compareChannelData cdd
    end
    clear compareChannelCount channelData
end
clear channelCount









