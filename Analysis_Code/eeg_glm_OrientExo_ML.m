% From eeg_comppac.m function in the PAC-Tools eeglab extension

function [pacval, beta] = eeg_glm_OrientExo_ML(Phase, Amp)
% General linear Model PAC method
     X = [cos(Phase), sin(Phase) ones(size(Phase))];                      % Building Matrix of regressors. Note : glmfit adds a column of 1s
    [beta,~, stats] = glmfit(X,Amp,'normal','constant','off');      % Fit the GLM
    
    % Remove warning
    w = warning('query','last');
    id = w.identifier;
    warning('off',id)
    
    pacval = 1- sum(stats.resid.^2)/sum((Amp-mean(Amp)).^2);  % 1-var(stats.resid)/var(amplitude); % Calculate the explained variance
end