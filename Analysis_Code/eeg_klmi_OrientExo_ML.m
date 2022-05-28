function [pacval, peakangle, bin_average, nbins] = eeg_klmi_OrientExo_ML(Phase, Amp, nbinskl)
% Compute the Kullback-Leibler modulation index
% 
%   For reference, see:
%   Tort, ABL, Komorowski, R, Eichenbaum, H, & Kopell, N (2010). Measuring
%   phase-amplitude coupling between neuronal oscillations of different
%   frequencies. J Neurophys, 104: 1195-1210.
% 
% 
% Code originally from PACTools eeglab extension: https://github.com/sccn/PACTools 
% Modified by SSS
    
    verbose = 0;
    modbin_flag = 1; %automatically adjust bin number when needed
    nbins = nbinskl;
    
    flag = true;
    while flag
        bin_size     = 2*pi/nbins;             % Compute bin size
        bin_average  = zeros(1,nbins);         % Initialize the histogramme

        for i = 1:nbins                         %Cycle through the bins   
        %   fill the bins with the mean amplitude for the phases that correspond to that bin
            bin_average(i) = mean(Amp(wrapTo2Pi(Phase) > (i-1)*bin_size & wrapTo2Pi(Phase) < i*bin_size));
        end
        clear i

        % Normalize the bins to obtain a probability distribution
        bin_average = bin_average/(sum(bin_average));

        % If some bins are empty, then the KL distance cannot be computed (log of a null value)
         if sum(bin_average>0)<nbins && modbin_flag
            oldnbins = nbins;
            nbins = ceil(nbins/2);
%             if verbose
%                 fprintf('Too many bins to compute Kullback Leibler MI. Reducing the number of bins from %d to %d \n', oldnbins, nbins);
%             end
         else
            flag = false;
         end
    end
    
    % Compute KL distance
    pacval = (sum(bin_average.*log(bin_average/(1/nbins))))/log(nbins);
    
    % Compute the peak phase
    [~,index] = max(bin_average);
    peakangle = index*bin_size - bin_size/2;
end