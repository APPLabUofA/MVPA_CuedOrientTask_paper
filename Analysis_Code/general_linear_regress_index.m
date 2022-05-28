%  Source: https://github.com/ChloeQin3104/laughing-dollop
% 
%   adapted from the code written for Kramer et al. (2008) by Tolga Ozkurt (2011)
%   to implement a General Linear Model (GLM) phase-amplitude coupling estimator

%   INPUTS:
%    a: amplitude, 
%    p: phase
%
%   OUTPUS:
%    
%    mod2d_raw1  = The two-dim GLM estimate
%    mod2d_raw2 =  The two-dim GLM estimate from regreession coefffients
%    mod2d_raw3 =  The two-dim spurious term removed robust GLM estimate 


 function [m_raw1, m_raw2, m_raw3] = general_linear_regress_index(Amp, Phase)

   % General Linear regression Model PAC estimator
   % a: amplitude, p: phase


numpoints = length(Amp);

X = [cos(Phase)' sin(Phase)' ones(numpoints,1)];
[beta_coef, ~, error_trms] = regress(Amp', X);
  
% standard GLM (Penny et al., 2008)
m_raw1 = sqrt( (sum(Amp.*Amp) - sum(error_trms.*error_trms)) / sum(Amp.*Amp) );
 
% equivalent to standard GLM shown by (Ozkurt and Schnitzler, 2011)
m_raw2 = sqrt ( numpoints * (beta_coef(1)^2 + beta_coef(2)^2 + beta_coef(3)^2) / sum(Amp.*Amp)  );
 
% robust GLM without the spurious term (beta_3) derived by (Ozkurt and Schnitzler, 2011)  
m_raw3 = 0.5 * sqrt ( (beta_coef(1)^2 + beta_coef(2)^2) / sum(Amp.*Amp)  );
 

mod2d_raw1 = abs(m_raw1);
mod2d_raw2 = abs(m_raw2);
mod2d_raw3 = abs(m_raw3);
