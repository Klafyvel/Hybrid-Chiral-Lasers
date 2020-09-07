%interface_chiral_same_handedness.m
%
%Calculate the matrix M such that $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_1 = \mathbf{M}\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_2$ where 1 is a chiral medium and 2 is a chiral medium characterized by $p$, $\epsilon_a$, $\epsilon_b$ 
%

function [M] = interface_chiral_same_handedness(p1, epsilona1, epsilonb1, psi1, p2, epsilona2, epsilonb2, psi2, z, k_0, lhm)
    if nargin < 11
        lhm = false;
    end

    [k1, kappa1, deltak1, n_bar1, deltan1] = cwt_convenient_variables(p1, epsilona1, epsilonb1, k_0,lhm);
    [k2, kappa2, deltak2, n_bar2, deltan2] = cwt_convenient_variables(p2, epsilona2, epsilonb2, k_0,lhm);

    M = interface_chiral_same_handedness_from_reduced_variables(p1,kappa1,n_bar1,psi1,p2,kappa2,n_bar2,psi2,z,k_0,lhm);
end
