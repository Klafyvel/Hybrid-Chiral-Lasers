% cwt_from_reduced_variables.m
% 
% Compute the propagation matrix from CW theory given reduced variables (useful for lasing plot)
% 
% The matrix is expressed in the $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}$ basis.
% 
% Parameters:
%     - p : parameter of the medium
%     - psi : phase of the medium
%     - L : length of the slab
%     - k : reduced wavelength
%     - kappa : Coupling constant
%     - deltak : Detuning
%     - left handedness (optional, default=false) : return the matrix for a left-handed medium.
function [M] = cwt_from_reduced_variables(p, psi, L, k, kappa, deltak, lhm)

if nargin < 7 
    lhm = false;
end

% kappaL = kappa*L;
% dkL = deltak*L;
% Dnu = L;

Delta = sqrt((kappa)^2-(deltak/2)^2);
% DeltaL        = sqrt((kappaL.*kappaL)-(dkL./2).*(dkL./2));

% M gives the propagation in the circular basis, i.e.
% [E⁺_L;E⁺_R;E⁻_L;E⁻_R]_{z=L} = M*[E⁺_L;E⁺_R;E⁻_L;E⁻_R]_{z=0}
if lhm
Q_plus  = i*exp(-i*(p*L+2*psi))*kappa/Delta*sinh(Delta*L);
Q_minus = -i*exp(i*(p*L+2*psi))*kappa/Delta*sinh(Delta*L);
P_plus  = exp(-i*p*L) * (cosh(Delta*L)+i*deltak*sinh(Delta*L)/(2*Delta));
P_minus = exp(i*p*L) * (cosh(Delta*L)-i*deltak*sinh(Delta*L)/(2*Delta));
M = [
P_plus  0           Q_plus  0
0       exp(i*k*L)  0       0
Q_minus 0           P_minus 0
0       0           0       exp(-i*k*L)
];
else
Q_plus  = i*exp(i*(p*L+2*psi))*kappa/Delta*sinh(Delta*L);
Q_minus = -i*exp(-i*(p*L+2*psi))*kappa/Delta*sinh(Delta*L);
P_plus = exp(i*p*L) * (cosh(Delta*L)+i*deltak*sinh(Delta*L)/(2*Delta));
P_minus = exp(-i*p*L) * (cosh(Delta*L)-i*deltak*sinh(Delta*L)/(2*Delta));
% P_plus = (cosh(DeltaL)+i.*(dkL./(2.*DeltaL)).*sinh(DeltaL)).*exp(i*p*Dnu);
% Q_plus         = (i.*kappaL./(DeltaL).*sinh(DeltaL)).*exp(i*p*Dnu);
% P_minus        = (cosh(DeltaL)-i.*(dkL./(2.*DeltaL)).*sinh(DeltaL)).*exp(-i*p*Dnu);
% Q_minus        = -(i.*kappaL./(DeltaL).*sinh(DeltaL)).*exp(-i*p*Dnu);
M = [
exp(i*k*L) 0       0           0
0          P_plus  0           Q_plus
0          0       exp(-i*k*L) 0
0          Q_minus 0           P_minus
];
end
end
