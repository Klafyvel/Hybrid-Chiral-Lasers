function [M] = interface_chiral_same_handedness_from_reduced_variables(p1,kappa1,n_bar1,psi1,p2,kappa2,n_bar2,psi2,z,k0,lhm)

    if nargin < 11
        lhm = false;
    end

    P1m = kappa1/k0 * exp(-2*i*(p1*z+psi1));
    P1p = kappa1/k0 * exp(2*i*(p1*z+psi1));
    P2m = kappa2/k0 * exp(-2*i*psi2);
    P2p = kappa2/k0 * exp(2*i*psi2);

    n_bar = (n_bar1+n_bar2)/2;
    deltan = n_bar2 - n_bar1;

    if ~lhm
        M = [
        4*n_bar2*n_bar-kappa2^2/k0^2 2*(n_bar*P2p-n_bar2*P1m) deltan*P2m 2*n_bar2*deltan+P1p*P2m-kappa2^2/k0^2;
        -deltan*P2p 4*n_bar2*n_bar-P1m*P2p 2*n_bar2*deltan 2*n_bar2*P1p-2*n_bar*P2p;
        deltan*P2p 2*n_bar2*deltan-kappa2^2/k0^2+P1m*P2p 4*n_bar2*n_bar-kappa2^2/k0^2 2*n_bar*P2p-2*n_bar2*P1p;
        2*n_bar2*deltan 2*n_bar2*P1m-2*n_bar*P2m -2*deltan*P2m 4*n_bar2*n_bar-P1p*P2m;
        ];
    else
        M = [
        4*n_bar2*n_bar-P1p*P2m -deltan*P2m 2*(n_bar2*P1m-n_bar*P2m) 2*n_bar2*deltan;
        2*(n_bar*P2p-n_bar2*P1p) 4*n_bar2*n_bar-kappa2^2/k0^2 2*n_bar2*deltan+P1m*P2p-kappa2^2/k0^2 deltan*P2p;
        2*(n_bar2*P1p-n_bar*P2p) 2*n_bar2*deltan 4*n_bar2*n_bar-P1m*P2p -deltan*P2p;
        2*n_bar2*deltan+P1p*P2m-kappa2^2/k0^2 deltan*P2m 2*(-n_bar2*P1m+n_bar*P2m) 4*n_bar*n_bar2-kappa2/k0^2;
        ];
    end
    M = 1/(4*n_bar2^2-kappa2^2/k0^2) * M; 
end
