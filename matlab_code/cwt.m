% cwt.m
% 
% Compute the M matrix (propagation matrix) using the Coupled Wave Theory for a chiral media defined by the parameters $p$, $\psi$, and $\epsilon_{a,b,c}$, where the permittivity is :
% $$
% \mathbf{\epsilon} = \mathbf{R}\cdot\begin{pmatrix}
% \epsilon_a  & 0          & 0 \\
% 0           & \epsilon_b & 0 \\
% 0           & 0          & \epsilon_c
% \end{pmatrix}\cdot\mathbf{R}^{-1}
% $$
% 
% with $\mathbf{R}$ being the matrix for a rotation of $pz+\psi$ around $z$ axis.
% The matrix is expressed in the $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}$ basis.
% 
% Parameters:
%     - p : parameter of the medium
%     - epsilona : $\epsilon_a$
%     - epsilonb : $\epsilon_b$
%     - psi : phase of the medium
%     - d : length of the slab
%     - k0 : reduced wavelength of the wave
%     - left handedness (optional, default=false) : return the matrix for a left-handed medium.
function [M] = cwt(p, epsilona, epsilonb, psi, d, k_0, lhm)

if nargin < 7 
    lhm = false;
end

[k, kappa, deltak] = cwt_convenient_variables(p, epsilona, epsilonb, k_0, lhm);

M = cwt_from_reduced_variables(p, psi, d, k, kappa, deltak, lhm);
end
