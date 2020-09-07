% hybrid_defect_tuning.m
%
% This file is used to simulate an hybrid chiral cavity and tune it.
%
% The cavity looks like this
%
%        interface 2         interface LR           interface RR       interface RL         interface 2
% isotropic  | LHCM               | RHCM               | RHCM               | LHCM               | isotropic 
% medium     | n_bar_L, delta_n_L | n_bar_R, delta_n_R | n_bar_R, delta_n_R | n_bar_L, delta_n_L | medium
% n_1        | p_L                | p_R                | psi, p_R           | psi, p_L           | n_2
%             <------ L_L -------> <-------- L/2------> <------- L/2-------> <------- L_L ------> 
close all;
clear all;

%% Data
Lp = 300e-9             % Spatial period of chiral media (same in LH and RH)

n_bar_L = 1.7690+0.02i;       % Average refractive index for LH chiral media (no gain)
kappaL_L = 4;           % Coupling constant * length for LH chiral media

% L_L = 50*Lp;            % Length of the DBRs
% L = 100*Lp;              % Length of RH medium
L_L = 40*Lp;            % Length of the DBRs
L = 20*Lp;              % Length of RH medium

re_n_bar_R = 1.7690;    % Real part of average refractive index in RH medium
kappaL_R = 4;           % Coupling constant * length for RH chiral media
max_detuning = 6;     % Maximum detuning in the RH medium Re(\delta k L / 2)
max_gain = 3;           % Maximum gain in RH medium -Im(\delta k L / 2)
psis = linspace(-pi, pi, 200); % Rotation of the RH medium

n1 = 1;%real(n_bar_L);  % Refractive index of isotropic medium 1
n2 = n1;                % Refractive index of isotropic medium 2

step_gain = 0.02;       % Step for the range of gains explored
step_detuning = 0.02;   % Step for the range of detunings explored

peak_radius = 0.2;      % Expected radius of lasing peaks in the lasing plot. 
                        % Reduce if you think some peaks are missing, increase 
                        % if you think there are false positives.

nb_col_mode = 1;        % Number of columns for the plot of modes.
width = 600;            % Width of the mode plots
%% Useful intermediary variables (Automatically set from data)

p = 2*pi/Lp;            % Periodicity of the chiral media
p_R = p;                % Periodicity for RH medium
p_L = -p;               % Periodicity for LH media
gain_range = -max_gain*0.1:step_gain:max_gain; % Gain range explored (-Im(\delta k \times L / 2))
detuning_range = 0:step_detuning:max_detuning; % Detuning range explored (Re(\delta k \times L / 2))

[X,Y] = meshgrid(detuning_range,gain_range); % X and Y meshes for detuning and gain in RHM

deltak_R = 2.*X ./ L - 2.*i.*Y./L; % Complex detuning in RHM
k_R = deltak_R ./ 2 + p_R;% Wavevector in RH medium

k_0 = real(k_R) ./ re_n_bar_R; % Wavevector in vacuum

n_bar_R = k_R ./ k_0;   % Refractive index of the RH medium, including gain
kappa_R = kappaL_R/L;   % Coupling constant of RH medium

kappa_L = kappaL_L/L_L; % Coupling constant in LH media
k_L = k_0 .* n_bar_L;   % Wavevector in LH medium
deltak_L = 2.*(k_L+p_L);  % Detuning in LH medium 

% permittivity in LH media, on each axes
epsilon_a_L = n_bar_L.^2 + (2 * kappa_L .* n_bar_L) ./ k_0;
epsilon_b_L = n_bar_L.^2 - (2 * kappa_L .* n_bar_L) ./ k_0;
% permittivity in RH medium, on each axes
epsilon_a_R = n_bar_R.^2 + (2 * kappa_R .* n_bar_R) ./ k_0;
epsilon_b_R = n_bar_R.^2 - (2 * kappa_R .* n_bar_R) ./ k_0;

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

%% Compute the invert determinant for the whole medium

% In order to improve speed, pre-allocate the array
lasing_coeff = zeros(length(gain_range), length(detuning_range));
% The cavity matrices will be stored here
cavities = zeros(4,4,length(gain_range), length(detuning_range));

row{length(psis)} = [];
col{length(psis)} = [];
lcp{length(psis)} = [];
epsilon{length(psis)} = [];

disp('Starting simulation');
for m=1:length(psis)
    psi = psis(m);
    for j=1:length(detuning_range)
        for l=1:length(gain_range)
            % Transfer matrix inside of the medium
            M_L1 = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, L_L, k_0(l,j));
            M_L2 = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), psi, L_L, k_0(l,j));
            M_R1 = oseen(p_R, epsilon_a_R(l,j), epsilon_b_R(l,j), 0, L/2, k_0(l,j));
            M_R2 = oseen(p_R, epsilon_a_R(l,j), epsilon_b_R(l,j), psi, L/2, k_0(l,j));
            % Transfer matrix, accounting for the interfaces, in circular basis.
            M_cavity = T^-1*M_int_2^-1*M_L2*M_R2*M_R1*M_L1*M_int_1*T;
            lasing_coeff(l,j) = abs(inv(det(M_cavity(3:4,3:4))));
            cavities(:,:,l,j) = M_cavity;
        end
    end
    % Finding lasing loci by looking for peaks in the log of the lasing coefficient
    coef_log = log(lasing_coeff);
    [row_m, col_m] = find_max(squeeze(coef_log(:,:)), peak_radius/step_detuning, peak_radius/step_gain);
    row{m} = row_m;
    col{m} = col_m;
    % Identify modes characteristics
    [R_m, theta_m, epsilon_m, lcp_m, modes_m] = find_modes(row_m,col_m,cavities);
    lcp{m} = lcp_m;
    epsilon{m} = epsilon_m;
    fprintf('Done:\t%d/%d\r', m,length(psis));
end
%%

row = horzcat(row{:});
col = horzcat(col{:});
lcp = horzcat(lcp{:});
epsilon = horzcat(epsilon{:});
%%
left_handed_modes = deltak_R(sub2ind(size(deltak_R), row(lcp).', col(lcp).'));
right_handed_modes = deltak_R(sub2ind(size(deltak_R), row(~lcp).', col(~lcp).'));

left_handed_gains = -imag(left_handed_modes*L/2);
left_handed_detunings = real(left_handed_modes*L/2);
right_handed_gains = -imag(right_handed_modes*L/2);
right_handed_detunings = real(right_handed_modes*L/2);


%%

min_epsilon = min(epsilon)
max_epsilon = max(epsilon)

%% Saving workspace

save('hybrid_defect_tuning')

%% Plotting

figure
scatter(right_handed_detunings, right_handed_gains, 'b.')
hold on;
scatter(left_handed_detunings, left_handed_gains, 'r.')
hold off;
legend('RCP', 'LCP', 'interpreter', 'latex')
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')

savefig(gcf, 'hybrid_defect_cavity', 'tuning', 'eps', 'epsc')

%% Plot purity

figure
histogram(epsilon*180/pi)
xlabel('$\epsilon (^\circ)$', 'interpreter', 'latex')
savefig(gcf, 'hybrid_defect_cavity', 'purity', 'eps', 'epsc')
figure
histogram(epsilon(lcp)*180/pi)
xlabel('$\epsilon (^\circ)$', 'interpreter', 'latex')
savefig(gcf, 'hybrid_defect_cavity', 'purity_lh', 'eps', 'epsc')
figure
histogram(epsilon(~lcp)*180/pi)
xlabel('$\epsilon (^\circ)$', 'interpreter', 'latex')
savefig(gcf, 'hybrid_defect_cavity', 'purity_rh', 'eps', 'epsc')

length(epsilon)
length(epsilon(lcp))
length(epsilon(~lcp))

%% Plot purity another way

figure
scatter(left_handed_detunings, left_handed_gains, 10, epsilon(lcp).'*180/pi, 'o', 'filled')
hold on
scatter(right_handed_detunings, right_handed_gains, 10, epsilon(~lcp).'*180/pi, 'o', 'filled')
hold off
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
zlabel('$\epsilon$', 'interpreter', 'latex')
colorbar()
savefig(gcf, 'hybrid_defect_cavity', 'purity_3D', 'eps', 'epsc')


%% Plot ellipses
[x, y] = ellipse(min_epsilon, 0, 0, 0);
plot(x,y, 'LineWidth', 2)
hold on
[x, y] = ellipse(max_epsilon, 0, 0, 0);
plot(x,y, '--', 'LineWidth', 2)
legend(sprintf('$\\epsilon=%.3g^{\\circ}$', min_epsilon*180/pi), sprintf('$\\epsilon=%.3g^{\\circ}$', max_epsilon*180/pi), 'interpreter', 'latex')
daspect([1 1 1])
hold off
savefig(gcf, 'hybrid_defect_cavity', 'ellipses', 'eps', 'epsc')
