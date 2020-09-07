% interface_chiral_to_isotrope_from_reduced_variables.m
%
% Calculate the matrix M such that $\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_1 = \mathbf{M}\begin{bmatrix}E^+_L\\E^+_R\\E^-_L\\E^-_R\end{bmatrix}_2$ where 1 is an isotrope medium of refractive index $n_1$ and 2 is a chiral medium characterized by $p$, $\epsilon_a$, $\epsilon_b$ 
%

function [M] = interface_chiral_to_isotrope_from_reduced_variables(p, psi, z, k, kappa, n1, k_0, lhm)

    if nargin < 8
        lhm = false;
    end

    n_bar = k / k_0;

    P_minus = kappa/(k_0*n1) * exp(-2*i*(p*z+psi));%deltan/(2*n1) * exp(-2*i*(p*z+psi));
    P_plus = kappa/(k_0*n1) * exp(2*i*(p*z+psi));%deltan/(2*n1) * exp(2*i*(p*z+psi));
    
    if lhm
        M = 1/2 * [
            1+n_bar/n1 0          P_minus    1-n_bar/n1;
            -P_plus    1+n_bar/n1 1-n_bar/n1 0;
            P_plus     1-n_bar/n1 1+n_bar/n1 0;
            1-n_bar/n1 0          -P_minus   1+n_bar/n1;
        ];
    else
        M = 1/2 * [
            1+n_bar/n1 -P_minus    0           1-n_bar/n1;
            0           1+n_bar/n1 1-n_bar/n1  P_plus;
            0           1-n_bar/n1 1+n_bar/n1 -P_plus;
            1-n_bar/n1  P_minus    0           1+n_bar/n1;
        ];
    end
end
