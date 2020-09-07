% oseen.m
% 
% Compute the M matrix (propagation matrix) using the Oseen transformation for a chiral media defined by the parameters $p$, $\psi$, and $\epsilon_{a,b,c}$, where the permittivity is :
% $$
% \mathbf{\epsilon} = \mathbf{R}\cdot\begin{pmatrix}
% \epsilon_a  & 0          & 0 \\
% 0           & \epsilon_b & 0 \\
% 0           & 0          & \epsilon_c
% \end{pmatrix}\cdot\mathbf{R}^{-1}
% $$
% 
% with $\mathbf{R}$ being the matrix for a rotation of $pz+\psi$ around $z$ axis.
% 
% The matrix is expressed in the $\begin{bmatrix}\mathbf{E}_{\perp}\\ \mathbf{H}_{\perp}\end{bmatrix}$ basis.
% 
% Parameters:
%     - p : parameter of the medium
%     - epsilona : $\epsilon_a$
%     - epsilonb : $\epsilon_b$
%     - psi : $\psi$
%     - d : length of the slab
%     - k0 : reduced wavelength
function [M] = oseen(p, epsilona, epsilonb, psi, d, k_0)
G_tilde = [
0            -i*p           0    k_0
i*p           0            -k_0  0 
0            -k_0*epsilonb  0   -i*p
k_0*epsilona  0             i*p  0
];
transformation = [rotation(psi+p*d) zeros(2, 2); zeros(2, 2) rotation(psi+p*d)];
transformation_inv = [rotation(-psi) zeros(2, 2); zeros(2, 2) rotation(-psi)];
M = transformation * expm(i*G_tilde*d) * transformation_inv;
end
