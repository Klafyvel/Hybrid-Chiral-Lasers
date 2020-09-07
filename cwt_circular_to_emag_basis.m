% cwt_circular_to_emag_basis.m
% 
% Compute the matrix to change a field written in the $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}$ basis into the $\begin{bmatrix}\mathbf{E}_{\perp}\\ \mathbf{H}'_{\perp}\end{bmatrix}$ basis.
% 
% For details on the parameters, see `cwt` doc.
% 
% Parameters:
%     - p : parameter of the medium
%     - epsilona : $\epsilon_a$
%     - epsilonb : $\epsilon_b$
%     - psi : $\psi$
%     - z : position of the calculation
%     - omega : reduced frequency of the wave

function [M] = cwt_circular_to_emag_basis(p, epsilona, epsilonb, psi, z, k0,lhm)
    if nargin < 7
        lhm = false;
    end
    
    [k, kappa, deltak, n_bar, deltan] = cwt_convenient_variables(p, epsilona, epsilonb, k0,lhm);
    
    M = cwt_circular_to_emag_basis_from_reduced_variables(p,psi,kappa,n_bar,z,k0,lhm);
end
