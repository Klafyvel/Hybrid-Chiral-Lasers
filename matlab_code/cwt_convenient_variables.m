% cwt_convenient_variables.m
% 
% Calculate the useful variables needed for CWT approach : $k$, $\kappa$ and $\delta k$.
% 
% Parameters:
%     - p : parameter of the medium
%     - epsilona : $\epsilon_a$
%     - epsilonb : $\epsilon_b$
%     - k_0 : reduced frequency of the wave
%     - lhm : if true the medium is left-handed

function [k, kappa, deltak, n_bar, deltan] = cwt_convenient_variables(p, epsilona, epsilonb, k_0, lhm)
if nargin < 5 
    lhm = false;
end
lambda_0 = 2 * pi / k_0  ;
n_tilde  = sqrt(epsilona)      ;
n_b      = sqrt(epsilonb)      ;
n_bar    = sqrt((epsilona+epsilonb)/2); %(n_b + n_tilde) / 2 ;
deltan   = (n_tilde - n_b)     ;
k        = n_bar * k_0         ;
delta_epsilon = (epsilona - epsilonb) / 2;
kappa    = (pi/lambda_0)*(delta_epsilon/n_bar); %pi*deltan/lambda_0  ;
if lhm
deltak   = 2*(k+p)             ;
else
deltak   = 2*(k-p)             ;
end
end
