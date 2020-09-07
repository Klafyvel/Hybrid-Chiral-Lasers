% Convert [E_L^+ E_R^+ E_L^- E_R^-]^T to [E_x E_y H_x H_y]^T in chiral medium

function [M] = cwt_circular_to_emag_basis_from_reduced_variables(p,psi,kappa,n_bar,z,k_0,lhm)
    if nargin < 7
        lhm = false;
    end
    P_plus = kappa/k_0*exp(i*2*(p*z+psi));
    P_minus = kappa/k_0*exp(-i*2*(p*z+psi));

    if ~lhm
        M = 1/sqrt(2) * [
        1        1                  1       1;
        i       -i                 -i       i;
        -i*n_bar i*(n_bar+P_minus) -i*n_bar i*(n_bar+P_plus);
        n_bar    n_bar-P_minus     -n_bar  -(n_bar-P_plus);
        ];
    else
        M = 1/sqrt(2) * [
        1                 1        1                 1;
        i                -i       -i                 i;
        -i*(n_bar+P_plus) i*n_bar -i*(n_bar+P_minus) i*n_bar;
        n_bar-P_plus      n_bar   -n_bar+P_minus    -n_bar;
        ];
    end
end
