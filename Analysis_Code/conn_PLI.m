%% phase_lag_index 
% This function computes the phase lag index (PLI) between the signals of 
% different channels or ROIs in the input matrix.
%
% PLI = phase_lag_index(sig)
%
% INPUT:
%   sig is the input matrix
%
% OUTPUT:
%   PLI is the matrix of the phase lag index (numTrials, numChannels, numChannels)
%
% Info:
%           Computes the phase lag index proposed by Stam et al., 2007 for 
%           each possible pair of EEG channels and for every band as well.
% 
%--------------------------------------------------------------------------
% NOTE:
% In order to extract the pli between channels 17 and 20, use pli(:,17,20) 
% and NOT pli(:,20,17). The smaller channel number is to be used first.
%--------------------------------------------------------------------------
%
% Mathematical background:
%           According to Stam et al., 2007 PLI is defined in a complete 
%           mathematical formula as:
%
%                   PLI = mean(|sign([sin(Df(t))])|)        (1)
%
%           where:
%
%                   Df(t) = phase1(t)-phase2(t)), is the phase difference
%                   of the two signals at time t, 
%
%                   phase1(t) is the phase of the 1st signal and is equal
%                   to arctan(x1_H(t)/x1(t)) where x1_H(t) is the Hilbert
%                   transformed version of signal x1(t),
%
%                   phase2(t) is the phase of the 2nd signal and is equal
%                   to arctan(x2_H(t)/x2(t)) where x2_H(t) is the Hilbert
%                   transformed version of signal x2(t).
%
% Fundamental basis:
%               The PLI ranges between 0 and 1. 
%               A PLI of zero indicates either no coupling or coupling with 
%               a phase difference centered around 0 mod p. A PLI of 1 
%               indicates perfect phase locking at a value of Df different 
%               from 0 mod p. The stronger this nonzero phase locking is, 
%               the larger PLI will be.
%
% Reference(s): Stam C. J., Nolte G., and Daffertshofer A (2007). Phase lag 
% index: assessment of functional connectivity from multi channel EEG and 
% MEG with diminished bias from common sources. Hum. Brain Mapp. 28(3),
% 1178-1193.            
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pli = conn_PLI(eegData) 

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

pli = NaN(numTrials, numChannels, numChannels); %pre-allocate
for channelCount = 1:numChannels-1
    channelData = squeeze(complex_sig(channelCount, :, :));
    for compareChannelCount = channelCount+1:numChannels
        compareChannelData = squeeze(complex_sig(compareChannelCount, :, :));
        
        % same equation as below
%         pli(:, channelCount, compareChannelCount) = ...
%             abs(mean(sign(angle(channelData./compareChannelData)),1));   
        
        % cross-spectral density
        cdd = channelData .* conj(compareChannelData);
        
        % phase-lag index
        pli(:, channelCount, compareChannelCount) = abs(mean(sign(imag(cdd)),1));
        
        clear compareChannelData cdd
    end
    clear compareChannelCount channelData
end
clear channelCount




