% reflection_cwt.m
% 
% Calculate the amplitude reflection matrix for a matrix M linking $\begin{bmatrix}E_L^+(d) \\E_R^+(d) \\E_L^-(d) \\E_R^-(d) \\ \end{bmatrix}$ and $\begin{bmatrix}E_L^+(0) \\E_R^+(0) \\E_L^-(0) \\E_R^-(0) \\\end{bmatrix}$.
% 

function [R] = reflection_cwt(M)
R = -inv(M(3:4,3:4))*M(3:4,1:2);
end
