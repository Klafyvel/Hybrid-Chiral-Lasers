% reflection_oseen.m
% 
function [T] = reflection_oseen(M_oseen, n1, n2)

% N1 = [
% 1 0 1 0;
% 0 1 0 1;
% 0 -n1 0 n1;
% n1 0 -n1 0;
% ];
% 
% N2 = [
% 1 0 1 0;
% 0 1 0 1;
% 0 -n2 0 n2;
% n2 0 -n2 0;
% ];

M_int_1 = [
    1 0 1 0;
    0 1 0 1;
    0 -n1 0 n1;
    n1 0 -n1 0;
];
M_int_2 = [
    1 0 1 0;
    0 1 0 1;
    0 -n2 0 n2;
    n2 0 -n2 0;
];
% Transfer matrix from circular to cartesian
T = 1/sqrt(2)*[
    1 1 0 0;
    i -i 0 0;
    0 0 1 1;
    0 0 -i i;
];


M = T^-1*M_int_2^-1*M_oseen*M_int_1*T;


T = M(1:2,1:2) - M(1:2,3:4)*M(3:4,3:4)^-1*M(3:4,1:2);
end
