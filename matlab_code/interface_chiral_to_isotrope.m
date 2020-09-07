% interface_chiral_to_isotrope.m
%
% Calculate the matrix M such that $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_1 = \mathbf{M}\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_2$ where 1 is an isotrope medium of refractive index $n_1$ and 2 is a chiral medium characterized by $p$, $\epsilon_a$, $\epsilon_b$ 
%

function [M] = interface_chiral_to_isotrope(p, epsilona, epsilonb, psi, z, k0, n1,lhm)
    if nargin < 8
        lhm = false;
    end
    
    [k, kappa, deltak, n_bar, deltan] = cwt_convenient_variables(p, epsilona, epsilonb, k0, lhm);

    M = interface_chiral_to_isotrope_from_reduced_variables(p,psi,z,k,kappa,n1,k0,lhm);

end
