% partial_inverse.m
% 
% computes the partial inverse of a matrix with pivot point at the center. See here : https://en.wikipedia.org/wiki/Partial_inverse_of_a_matrix

function [pM] = partial_inverse(M)
[m,n] = size(M);
M11 = M(1:m/2, 1:n/2);
M12 = M(1:m/2, n/2+1:end);
M21 = M(m/2+1:end, 1:n/2);
M22 = M(m/2+1:end, n/2+1:end);
M11_inv = M11^-1;
pM = [
M11_inv     -M11_inv*M12;
M21*M11_inv M22-M21*M11_inv*M12;
];
end
