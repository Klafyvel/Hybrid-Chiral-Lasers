% defect_tuning.m
%
% This file is used to simulate an defect chiral cavity and tune it.
%
% The cavity looks like this
%
%        interface 1       interface LL           interface 2
% isotropic  | LHCM             | LHCM               | isotropic 
% medium     | n_bar, delta_n   | n_bar_L, delta_n_L | medium
% n_1        | p                | psi, p             | n_2
%             <------ L -------> <-------- L -------> 
close all;
clear all;

%% Data

Lp = 300e-9             % Spatial period of chiral media (same in LH and RH)

% L = 100*Lp;              % Length of each slab of medium
L = 20*Lp;              % Length of each slab of medium

re_n_bar = 1.7690;      % Real part of average refractive index in RH medium
kappaL = 4;             % Coupling constant * length for RH chiral media
max_detuning = 12;       % Maximum detuning in the RH medium Re(\delta k L / 2)
max_gain = 2.5;           % Maximum gain in RH medium -Im(\delta k L / 2)
psis = linspace(-pi, pi, 200);                % Rotation of the RH medium

n1 = 1;                 % Refractive index of isotropic medium 1
n2 = n1;                % Refractive index of isotropic medium 2

step_gain = 0.02;       % Step for the range of gains explored
step_detuning = 0.02;   % Step for the range of detunings explored

peak_radius = 0.2;      % Expected radius of lasing peaks in the lasing plot. 
                        % Reduce if you think some peaks are missing, increase 
                        % if you think there are false positives.

nb_col_mode = 2;        % Number of columns for the plot of modes.

lhm = false;            % True if working with left-handed media.
width = 600;            % Width of the mode analysis.

%% Useful intermediary variables (Automatically set from data)

if lhm
    p = -2*pi/Lp;           % Periodicity of the chiral media
else
    p = 2*pi/Lp;
end
gain_range = -max_gain*0.1:step_gain:max_gain; % Gain range explored (-Im(\delta k \times L / 2))
detuning_range = 0:step_detuning:max_detuning; % Detuning range explored (Re(\delta k \times L / 2))

[X,Y] = meshgrid(detuning_range,gain_range); % X and Y meshes for detuning and gain in RHM

deltak = 2.*X ./ L - 2.*i.*Y./L; % Complex detuning in medium

if lhm
    k = deltak ./ 2 - p;    % Wavevector in LH medium
else
    k = deltak ./ 2 + p;    % Wavevector in RH medium
end

k_0 = real(k) ./ re_n_bar; % Wavevector in vacuum

n_bar = k ./ k_0;       % Refractive index of the RH medium, including gain
kappa = kappaL/L;       % Coupling constant of RH medium

% permittivity in LH media, on each axes
epsilon_a = n_bar.^2 + (2 * kappa .* n_bar) ./ k_0;
epsilon_b = n_bar.^2 - (2 * kappa .* n_bar) ./ k_0;

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

% In order to improve speed, pre-allocate the arrays
lasing_coeff = zeros(length(gain_range), length(detuning_range));
cavities = zeros(4,4,length(gain_range), length(detuning_range));

row{length(psis)} = [];
col{length(psis)} = [];
lcp{length(psis)} = [];

disp('Starting simulation');
for m=1:length(psis)
    psi = psis(m);
    for j=1:length(detuning_range)
        for l=1:length(gain_range)
            % Transfer matrix inside of the medium
            M_1 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), 0, L/2, k_0(l,j));
            M_2 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), psi, L/2, k_0(l,j));
            % Transfer matrix, accounting for the interfaces, in circular basis.
            M_cavity = T^-1*M_int_2^-1*M_2*M_1*M_int_1*T;
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
    disp(sprintf('Done:\t%d/%d', m,length(psis)));
end


row = horzcat(row{:});
col = horzcat(col{:});
lcp = horzcat(lcp{:});

left_handed_modes = deltak(sub2ind(size(deltak), row(lcp).', col(lcp).'));
right_handed_modes = deltak(sub2ind(size(deltak), row(~lcp).', col(~lcp).'));

left_handed_gains = -imag(left_handed_modes*L/2);
left_handed_detunings = real(left_handed_modes*L/2);
right_handed_gains = -imag(right_handed_modes*L/2);
right_handed_detunings = real(right_handed_modes*L/2);

%% Plotting

figure
scatter(right_handed_detunings, right_handed_gains, 'b.')
hold on;
scatter(left_handed_detunings, left_handed_gains, 'r.')
hold off;
legend('RCP', 'LCP', 'interpreter', 'latex')
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')

savefig(gcf, 'defect_cavity', 'tuning', 'eps', 'epsc')
