% find_modes.m
%
% Calculate the output modes for a given set of loci

function [R, theta, epsilon, lcp, modes] = find_modes(row,col,cavities)
% The modes at z=0 will be stored here
modes = zeros(2, length(row));

for u=1:length(row)
    l = row(u);
    j = col(u);

    M = cavities(:,:,l,j);

    [V,D] = eig(M(3:4,3:4)); % V is the matrix of eignevectors and D the diagonalized matrix

    % The corresponding eigen vector has an eigenvalue of nearly 0
    [m,min_index] = min(diag(D));
    
    % Use boundary condition to determine the full field at interface
    back_field = V(1:2,min_index);
    E0_L_minus = back_field(1);
    E0_R_minus = back_field(2);
    
    modes(:, u) = [E0_L_minus;E0_R_minus];   % The output mode in medium 1
end

% calculate the parameters for the ellipse
E_L = modes(1,:);
E_R = modes(2,:);
E_x = E_L + E_R;
E_y = i.*E_L - i.*E_R;

R = atan2(abs(E_x), abs(E_y));
delta = mod(angle(E_y), 2*pi) - mod(angle(E_x), 2*pi);
theta = 1/2 .* atan(tan(2.*R).*cos(delta));
epsilon = 1/2 .* asin(sin(2.*R).*abs(sin(delta)));
lcp = abs(E_L) > abs(E_R);
end
