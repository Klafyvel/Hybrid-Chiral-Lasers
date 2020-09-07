% interface_left_to_right_from_reduced_variables.m
% 
% This file computes the matrix associated to the inerface from a left-handed chiral medium to a right-handed chiral medium.

function [M] = interface_left_to_right_from_reduced_variables(pl, psil, kl, kappal, zl, pr, psir, kr, kappar, zr, k0)

n_barr = kr / k0;
n_barl = kl / k0;
deltanr = 2*kappal/k0;

P_minusl = kappal/k0 * exp(-2*i*(pl*zl+psil));
P_minusr = kappar/k0 * exp(-2*i*(pr*zr+psir));
P_plusl = kappal/k0 * exp(2*i*(pl*zl+psil));
P_plusr = kappar/k0 * exp(2*i*(pr*zr+psir));

n_bar = (n_barl+n_barr)/2;
deltan = n_barl - n_barr;

M = 1/(4*n_barr^2-(kappar/k0)^2)*[
    4*n_barr*n_bar-P_plusl*P_minusr-(kappar/k0)^2 2*n_bar*P_minusr             2*n_barr*P_minusl-deltan*P_minusr            -2*n_barr*deltan-(kappar/k0)^2;
    deltan*P_plusr-2*n_barr*P_plusl             4*n_barr*n_bar              -2*n_barr*deltan+P_minusl*P_plusr            -2*n_bar*P_plusr;
    2*n_barr*P_plusl-deltan*P_plusr            -2*n_barr*deltan-(kappar/k0)^2  4*n_barr*n_bar-P_minusl*P_plusr-(kappar/k0)^2 2*n_bar*P_plusr;
    -2*n_barr*deltan+P_minusr*P_plusl          -2*n_bar*P_minusr             deltan*P_minusr-2*n_barr*P_minusl           4*n_barr*n_bar;

];

end
