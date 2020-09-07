% interface_left_to_right.m
% 
% This file computes the matrix associated to the inerface from a left-handed chiral medium to a right-handed chiral medium.

function [M] = interface_left_to_right(pl, epsilonal, epsilonbl, psil, zl, pr, epsilonar, epsilonbr, psir, zr, k0)
    [kl, kappal, deltakl, n_barl, deltanl] = cwt_convenient_variables(pl, epsilonal, epsilonbl, k0);
    [kr, kappar, deltakr, n_barr, deltanr, k_0r] = cwt_convenient_variables(pr, epsilonar, epsilonbr, k0);

    M = interface_left_to_right_from_reduced_variables(pl, psil, kl, kappal, zl, pr, psir, kr, kappar, zr, k0);
end
