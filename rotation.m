% rotation.m
% 
% Compute the matrix $$\begin{pmatrix}\cos\ksi & s\sin\kis\\ \sin\kis & \cos\ksi \end{patrix}$$.
% 
% Parameters:
%     - ksi : rotation angle
% 
function [R] = rotation(ksi)
R = [
    cos(ksi) -sin(ksi)
    sin(ksi)  cos(ksi)
];
end
