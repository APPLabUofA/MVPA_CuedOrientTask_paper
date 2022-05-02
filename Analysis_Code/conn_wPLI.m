%% Weighted Phase Lag Index
% 
% [wpli,dwpli] = conn_wPLI(eegData) 
% 
% Input(s):
%           eegData  - input filtered dataset (channel x trials x time)
% 
% Outputs:
%           wpli     - weighted phase lag index matrix (trials x channel x
%                      channel)
%           dwpli    - debiased weighted phase-lag index matrix (trials x 
%                      channel x channel)
%--------------------------------------------------------------------------
% NOTE:
% In order to extract the wpli between channels 17 and 20, use wpli(:,17,20) 
% and NOT wpli(:,20,17). The smaller channel number is to be used first.
%--------------------------------------------------------------------------
% 
% Computes the weighted phase lag index from a data-matrix containing a 
% cross-spectral density. It implements the method described in: Vinck M, 
% Oostenveld R, van Wingerden M, Battaglia F, Pennartz CM. An improved 
% index of phase-synchronization for electrophysiological data in the 
% presence of volume-conduction, noise and sample-size bias. Neuroimage. 
% 2011 Apr 15;55(4):1548-65.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wpli,dwpli] = conn_wPLI(eegData) 

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

wpli = NaN(numTrials, numChannels, numChannels); %pre-allocate
dwpli = NaN(numTrials, numChannels, numChannels); %pre-allocate
for channelCount = 1:numChannels-1
    channelData = squeeze(complex_sig(channelCount, :, :));
    for compareChannelCount = channelCount+1:numChannels
        compareChannelData = squeeze(complex_sig(compareChannelCount, :, :));
                 
        % cross-spectral density
        cdd = channelData .* conj(compareChannelData);
        
        % take imaginary part of signal only
        cdi = imag(cdd);
        
        % weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
        wpli(:,channelCount,compareChannelCount) = abs(mean(abs(cdi).*sign(cdi),1))./mean(abs(cdi),1); % estimator of E(Im(X))/E(|Im(X)|)
        
        % debiased weighted phase-lag index (shortcut, as implemented in fieldtrip)
        imagsum        = sum(cdi,1); % compute the sum;
        imagsumW       = sum(abs(cdi),1); % normalization of the WPLI
        debiasfactor   = sum(cdi.^2,1);
        dwpli(:,channelCount,compareChannelCount) = mean((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor),1); % do the pairwise thing in a handy way
        
        clear compareChannelData cdd cdi imagsum imagsumW debiasfactor
    end
    clear compareChannelCount channelData
end
clear channelCount complex_sig



















