% hybrid_cavity.m
%
% This file is used to simulate an hybrid chiral cavity.
%
% The cavity looks like this
%
%        interface 2         interface LR           interface RL         interface 2
% isotropic  | LHCM               | RHCM                 | LHCM               | isotropic 
% medium     | n_bar_L, delta_n_L | n_bar_R, delta_n_R   | n_bar_L, delta_n_L | medium
% n_1        | p_L                | psi_R, p_R           | p_L                | n_2
%             <------ L_L -------> <-------- L ---------> <------- L_L ------> 
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
psi = 0;                % Rotation of the RH medium

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

% Positions for scanning intensity in LH media
intensity_scan_pos_L = linspace(0,L_L,400);
% Positions for scanning intensity in RH medium
intensity_scan_pos_R = linspace(0,L,2000);

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

% In order to improve speed, pre-allocate the array
lasing_coeff_oseen = zeros(length(gain_range), length(detuning_range));
% The cavity matrices will be stored here
cavities_oseen = zeros(4,4,length(gain_range), length(detuning_range));

% Scalar case
lasing_coeff_scalar = zeros(length(gain_range), length(detuning_range));

for j=1:length(detuning_range)
    for l=1:length(gain_range)
        % Matrix for right handed medium
        M_RHM = cwt_from_reduced_variables(p_R, psi, L, k_R(l,j), kappa_R, deltak_R(l,j));
        % Matrix for the two LH media
        M_LHM = cwt_from_reduced_variables(p_L, 0, L_L, k_L(l,j), kappa_L, deltak_L(l,j), true);
        % Matrix for interface 1 (isotrope -> LH)
        M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p_L, 0, 0, k_L(l,j), kappa_L, n1, k_0(l,j), true)^-1;
        % Matrix for interface LR (left handed -> right handed). Note that in the LH referencial 
        % this is at position L_L and in the RH referencial this is at position 0.
        M_interface_LR = interface_left_to_right_from_reduced_variables(p_L, 0, k_L(l,j), kappa_L, L_L, p_R, psi, k_R(l,j), kappa_R, 0, k_0(l,j));
        % Matrix for interface RL (right handed -> left handed). Note that in the LH referencial 
        % this is at position 0 and in the RH referencial this is at position L.
        M_interface_RL = interface_left_to_right_from_reduced_variables(p_L, 0, k_L(l,j), kappa_L, 0, p_R, psi, k_R(l,j), kappa_R, L, k_0(l,j))^-1;
        % Matrix for interface 2 (LH -> isotrope)
        M_interface_2 = interface_chiral_to_isotrope_from_reduced_variables(p_L, 0, L_L, k_L(l,j), kappa_L, n2, k_0(l,j), true);

        M_cavity = M_interface_2*M_LHM*M_interface_RL*M_RHM*M_interface_LR*M_LHM*M_interface_1;
        lasing_coeff(l,j) = abs(inv(det(M_cavity(3:4,3:4))));
        cavities(:,:,l,j) = M_cavity;
        cavities_inside(:,:,l,j) = M_RHM;

        % Oseen method
        % Transfer matrix inside of the medium
        M_L = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, L_L, k_0(l,j));
        M_R = oseen(p_R, epsilon_a_R(l,j), epsilon_b_R(l,j), 0, L, k_0(l,j));
        % Transfer matrix, accounting for the interfaces, in circular basis.
        M_oseen = T^-1*M_int_2^-1*M_L*M_R*M_L*M_int_1*T;
        % Here I use the boundary condition that no field is coming from outside
        lasing_coeff_oseen(l,j) = abs(inv(det(M_oseen(3:4,3:4))));
        cavities_oseen(:,:,l,j) = M_oseen;

        % Scalar case

        % Consider the two "mirrors" are identical, and calculate the reflection coefficient associated with each.
        M_interface_1 = interface_chiral_to_isotrope(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, 0, k_0(l,j), n1, true)^-1;
        M_medium = cwt(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, L_L, k_0(l,j), true);
        M_mirror = inv(M_interface_LR*M_medium*M_interface_1^-1);
        R = reflection_cwt(M_mirror);
        rLL = R(1,1);
        
        % The lasing coefficient in the scalar case takes into account the two reflections on the mirrors
        % and the round-trip in the cavity.
        lasing_coeff_scalar(l,j) = abs(inv(1-rLL^2*exp(2*i*k_R(l,j)*L_L)));
    end
end

%% Finding lasing loci by looking for pikes in the log of the lasing coefficient
% The modes are identified by peaks on the logarithm of the invert of the 
% determinant. 
coef_log = log(lasing_coeff);
% row and col will store the rows and columns where lasing modes are found
[row, col] = find_max(coef_log, peak_radius/step_detuning, peak_radius/step_gain);

coef_log_oseen = log(lasing_coeff_oseen);
% row_oseen and col_oseen will store the rows and columns where lasing modes are 
% found with oseen transformation 
[row_oseen, col_oseen] = find_max(coef_log_oseen, peak_radius/step_detuning, peak_radius/step_gain);

% Lasing coefficients and lasing loci for scalar case.
coef_log_scalar = log(lasing_coeff_scalar);
[row_scalar, col_scalar] = find_max(coef_log_scalar, peak_radius/step_detuning, peak_radius/step_gain);

%% Calculate the field at z=0 for each mode 

% R, theta and epsilon define the ellipse, lcp the handedness (true if left-handed) and modes are the eigen vectors, i.e. the output modes.
[R, theta, epsilon, lcp, modes] = find_modes(row,col,cavities);

% Same but for the Oseen method.
[R_oseen, theta_oseen, epsilon_oseen, lcp_oseen, modes_oseen] = find_modes(row_oseen,col_oseen,cavities_oseen);

%% Calculate the intensity distribution

% Number of scan positions in LH media
nb_scan_pos_L = length(intensity_scan_pos_L);
% Number of scan positions in RH medium
nb_scan_pos_R = length(intensity_scan_pos_R);

% Number of scan positions in the cavity
nb_scan_pos = 2*nb_scan_pos_L + nb_scan_pos_R;

% the intensities will be stored here
intensities = zeros(nb_scan_pos, 4*length(row));

for u=1:length(row)
    l = row(u);
    j = col(u);

    m = modes(:,u);

    % Field at z=0⁻
    F_0 = [0;0;m(1);m(2)];

    % Transfer matrix at interface 1
    M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p_L, 0, 0, k_L(l,j), kappa_L, n1, k_0(l,j), true)^-1;
    % Transfer matrix at interface Left to right
    M_interface_LR = interface_left_to_right_from_reduced_variables(p_L, 0, k_L(l,j), kappa_L, L_L, p_R, psi, k_R(l,j), kappa_R, 0, k_0(l,j));
    % Transfer matrix at interface right to left
    M_interface_RL = interface_left_to_right_from_reduced_variables(p_L, 0, k_L(l,j), kappa_L, 0, p_R, psi, k_R(l,j), kappa_R, L, k_0(l,j))^-1;

    % Shortcuts to store the results of intensity scan at the right place
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan the first LH medium
    for v=1:nb_scan_pos_L
        % length of the partial slab
        len = intensity_scan_pos_L(v);

        % Corresponding matrix
        M_LHM = cwt_from_reduced_variables(p_L, 0, len, k_L(l,j), kappa_L, deltak_L(l,j), true);
        % Field at z
        F = M_LHM * M_interface_1 * F_0;

        % Intensity at z
        intensities(v,res_start:res_stop) = abs(F).^2;
    end
    % Shortcut to store the matrix corresponding to the first LH medium
    M = M_interface_LR * cwt_from_reduced_variables(p_L, 0, L_L, k_L(l,j), kappa_L, deltak_L(l,j), true) * M_interface_1;

    % Scan the RH medium
    for v=1:nb_scan_pos_R
        % length of the partial slab
        len = intensity_scan_pos_R(v);

        % Corresponding matrix
        M_RHM = cwt_from_reduced_variables(p_R, psi, len, k_R(l,j), kappa_R, deltak_R(l,j));
        % Field at z
        F = M_RHM*(M * F_0);
        % Intensity at z
        intensities(v+nb_scan_pos_L,res_start:res_stop) = abs(F).^2;
    end
    % Shortcut to store the matrix corresponding to the first LH medium and the RH medium
    M = M_interface_RL * cwt_from_reduced_variables(p_R, psi, L, k_R(l,j), kappa_R, deltak_R(l,j)) * M;

    % Scan the second LH medium
    for v=1:nb_scan_pos_L
        % length of the partial slab
        len = intensity_scan_pos_L(v);

        % Corresponding matrix
        M_LHM = cwt_from_reduced_variables(p_L, 0, len, k_L(l,j), kappa_L, deltak_L(l,j), true);
        % Field at z
        F = M_LHM*(M * F_0);
        % Intensity at z
        intensities(v+nb_scan_pos_L+nb_scan_pos_R,res_start:res_stop) = abs(F).^2;
    end
end

%% Scan through the cavity with the Oseen method to check SVEA

% Pre-allocate the array for storing intensities
intensities_oseen = zeros(nb_scan_pos, 4*length(row_oseen));

for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);

    m = modes_oseen(:,u);

    % Field at z=0⁻ (electromagnetic basis)
    F_0 = M_int_1*T*[0;0;m(1);m(2)];

    % Shorthand to access the array of intensities at the right position
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan through first left-handed medium
    for v=1:nb_scan_pos_L
        % position in the medium, -> length of the virtual slab 
        len = intensity_scan_pos_L(v);
        % Corresponding matrix
        M_L = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, len, k_0(l,j));
        % Intensity at z
        intensities_oseen(v,res_start:res_stop) = abs(M_L*F_0).^2;
    end
    % Shorthand for the first LH medium
    M_L1 = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, L_L, k_0(l,j));

    % Scan through RH medium
    for v=1:nb_scan_pos_R
        % Position in RH medium
        len = intensity_scan_pos_R(v);
        % Corresponding matrix
        M_R = oseen(p_R, epsilon_a_R(l,j), epsilon_b_R(l,j), 0, len, k_0(l,j));
        % Intensity at z
        intensities_oseen(v+nb_scan_pos_L,res_start:res_stop) = abs(M_R*M_L1*F_0).^2;
    end
    % Shorthand for LH+RH medium
    M_R1 = oseen(p_R, epsilon_a_R(l,j), epsilon_b_R(l,j), 0, L, k_0(l,j)) * M_L1;

    % Scan through second LH medium
    for v=1:nb_scan_pos_L
        % Position in LH medium
        len = intensity_scan_pos_L(v);
        % Corresponding matrix
        M_L = oseen(p_L, epsilon_a_L(l,j), epsilon_b_L(l,j), 0, len, k_0(l,j));
        % Intensity at z
        intensities_oseen(v+nb_scan_pos_L+nb_scan_pos_R,res_start:res_stop) = abs(M_L*M_R1*F_0).^2;
    end
end

%% %%%%%%%%% BEYOND THAT POINT, THERE ARE ONLY PLOTS. %%%%%%%%% %%

%% Plot the surface and contour of lasing coefficient

figure('Position', [0 0 1200 600])
subplot(121)
surf(X,Y,coef_log, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Developped CWT')
subplot(122)
contour(X,Y,coef_log, 30)
xlabel('Detuning')
ylabel('Gain')
title('Developped CWT')
for j=1:length(row)
    x = X(row(j), col(j));
    y = Y(row(j), col(j));
    s = sprintf('< %.3g %+.3g i', x, y);
    text(x, y, s)
end

savefig(gcf, 'hybrid_cavity', 'lasing', 'png')

% Same thing but in EPS for the report.
figure('visible','off')
surf(X,Y,coef_log, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'hybrid_cavity', 'surface', 'eps', 'epsc')

figure('visible','off')
contour(X,Y,coef_log, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row)
    x = X(row(j), col(j));
    y = Y(row(j), col(j));
    s = sprintf('$\\leftarrow %.3g %+.3g i$', x, y);
    text(x, y, s,'interpreter','latex')
end
savefig(gcf, 'hybrid_cavity', 'contour', 'eps', 'epsc')


%% Plot the surface and contour of lasing coefficient (oseen)

figure('Position', [0 0 1200 600])
subplot(121)
surf(X,Y,coef_log_oseen, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Oseen Transformation')
subplot(122)
contour(X,Y,coef_log_oseen, 30)
xlabel('Detuning')
ylabel('Gain')
title('Oseen Transformation')
for j=1:length(row_oseen)
    x = X(row_oseen(j), col_oseen(j));
    y = Y(row_oseen(j), col_oseen(j));
    s = sprintf('< %.3g %+.3g i', x, y);
    text(x,y,s)
end

savefig(gcf, 'hybrid_cavity', 'lasing_oseen', 'png')

% Same thing but in EPS for the report.
figure('visible','off')
surf(X,Y,coef_log_oseen, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'hybrid_cavity', 'surface_oseen', 'eps', 'epsc')

figure('visible','off')
contour(X,Y,coef_log, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row_oseen)
    x = X(row_oseen(j), col_oseen(j));
    y = Y(row_oseen(j), col_oseen(j));
    s = sprintf('$\\leftarrow %.3g %+.3g i$', x, y);
    text(x, y, s,'interpreter','latex')
end
savefig(gcf, 'hybrid_cavity', 'contour_oseen', 'eps', 'epsc')

%% Plot the surface and contour of lasing coefficient (scalar)

figure('Position', [0 0 1200 600])
subplot(121)
surf(X,Y,coef_log_scalar, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Scalar case')
subplot(122)
contour(X,Y,coef_log_scalar, 30)
xlabel('Detuning')
ylabel('Gain')
title('Scalar case')
for j=1:length(row_scalar)
    x = X(row_scalar(j), col_scalar(j));
    y = Y(row_scalar(j), col_scalar(j));
    s = sprintf('< %.3g %+.3g i', x, y);
    text(x,y,s)
end

savefig(gcf, 'hybrid_cavity', 'lasing_scalar', 'png')

% Same thing but in EPS for the report.
figure('visible','off')
surf(X,Y,coef_log_scalar, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'hybrid_cavity', 'lasing_scalar_surf', 'eps', 'epsc')

figure('visible','off')
contour(X,Y,coef_log_scalar, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row_scalar)
    x = X(row_scalar(j), col_scalar(j));
    y = Y(row_scalar(j), col_scalar(j));
    s = sprintf('$\\leftarrow %.3g %+.3g i$', x, y);
    text(x, y, s,'interpreter','latex')
end
savefig(gcf, 'hybrid_cavity', 'lasing_scalar_contour', 'eps', 'epsc')

%% Plot intensity distribution

n_rows = ceil(length(row) / nb_col_mode);
if mod(length(row), nb_col_mode) == 0
    n_rows = n_rows + 1;
end
figure('Name', 'Intensity distribution for found modes', 'Position', [0 0 width (n_rows+2)*250])
subplot(n_rows+2,nb_col_mode,[1 2*nb_col_mode])
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)))
title('Identification of output modes, developped CWT')
xlabel('Detuning')
ylabel('Gain')
daspect([1 1 1])

X_intensities = cat(2, intensity_scan_pos_L, intensity_scan_pos_R+L_L, intensity_scan_pos_L+L+L_L) ./ (2*L_L+L);

for u=1:length(row)
    l = row(u);
    j = col(u);
    name = sprintf('\\textbf{Mode %d}\n($\\lambda_0$ = %.3g nm, $\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$)', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows+2,nb_col_mode,u+2*nb_col_mode)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, res_start), '-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    h_fl = fill([0 L_L/(L+2*L_L) L_L/(L+2*L_L) 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fl = fill([(L+L_L)/(L+2*L_L) 1 1 (L+L_L)/(L+2*L_L) ], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fr = fill([L_L/(L+2*L_L) (L+L_L)/(L+2*L_L) (L+L_L)/(L+2*L_L) L_L/(L+2*L_L)], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'-','LineWidth',2,'Color', [247 216 79]/255);
    h2 = plot(X_intensities, intensities(:, res_start), 'r-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), 'b-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows+2,nb_col_mode,[u+1+2*nb_col_mode nb_col_mode*(n_rows+2)]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h2;h3;h4;h5;h_fl;h_fr],'I','I_L^+','I_R^+','I_L^-','I_R^-','LH medium', 'RH medium','NumColumns',3);
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'hybrid_cavity', 'intensity_distribution', 'png')

% Save to EPS for the report
figure('Name', 'Intensity distribution for found modes', 'visible', 'off')
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)))
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
daspect([1 1 1])
savefig(gcf, 'hybrid_cavity', 'modes_found', 'eps', 'epsc')

figure('Name', 'Intensity distribution for found modes', 'Position', [0 0 width n_rows*250], 'visible', 'off')
for u=1:length(row)
    l = row(u);
    j = col(u);
    name = sprintf('\\textbf{Mode %d}\n($\\lambda_0$ = %.3g nm, $\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$)', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows,nb_col_mode,u)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, res_start), '-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    h_fl = fill([0 L_L/(L+2*L_L) L_L/(L+2*L_L) 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fl = fill([(L+L_L)/(L+2*L_L) 1 1 (L+L_L)/(L+2*L_L) ], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fr = fill([L_L/(L+2*L_L) (L+L_L)/(L+2*L_L) (L+L_L)/(L+2*L_L) L_L/(L+2*L_L)], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h2 = plot(X_intensities, intensities(:, res_start), 'r-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), 'b-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'-','LineWidth',2,'Color', [247 216 79]/255);
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows,nb_col_mode,[u+1 nb_col_mode*n_rows]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h2;h3;h4;h5;h_fl;h_fr],'$I$','$I_L^+$','$I_R^+$','$I_L^-$','$I_R^-$','LH medium', 'RH medium','NumColumns',3, 'interpreter', 'latex');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'hybrid_cavity', 'intensity_distribution', 'eps', 'epsc')

%% Plot intensity distribution from the Oseen method

n_rows = ceil(length(row_oseen) / nb_col_mode);
if mod(length(row_oseen), nb_col_mode) == 0
    n_rows = n_rows + 1;
end
figure('Name', 'Field distribution for found modes', 'Position', [0 0 width (n_rows+2)*250])
subplot(n_rows+2,nb_col_mode,[1 2*nb_col_mode])
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
title('Identification of output modes, Oseen method')
xlabel('Detuning')
ylabel('Gain')
daspect([1 1 1])

X_intensities = cat(2, intensity_scan_pos_L, intensity_scan_pos_R+L_L, intensity_scan_pos_L+L+L_L) ./ (2*L_L+L);

for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);
    name = sprintf('\\textbf{Mode %d}\n($\\lambda_0$ = %.3g nm, $\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$)', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows+2,nb_col_mode,u+2*nb_col_mode)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, intensities_oseen(:, res_start)+intensities_oseen(:, res_start+1));
    hold on
    h_fl = fill([0 L_L/(L+2*L_L) L_L/(L+2*L_L) 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fl = fill([(L+L_L)/(L+2*L_L) 1 1 (L+L_L)/(L+2*L_L) ], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fr = fill([L_L/(L+2*L_L) (L+L_L)/(L+2*L_L) (L+L_L)/(L+2*L_L) L_L/(L+2*L_L)], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h1 = plot(X_intensities, intensities_oseen(:, res_start)+intensities_oseen(:, res_start+1));
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows+2,nb_col_mode,[u+1+2*nb_col_mode nb_col_mode*(n_rows+2)]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h_fl;h_fr],'I','LH medium', 'RH medium');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'hybrid_cavity', 'intensity_distribution_oseen', 'png')

% Save to EPS for the report
figure('Name', 'Intensity distribution for found modes', 'visible', 'off')
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
daspect([1 1 1])
savefig(gcf, 'hybrid_cavity', 'modes_found_oseen', 'eps', 'epsc')

figure('Name', 'Intensity distribution for found modes', 'Position', [0 0 width n_rows*250], 'visible', 'off')
for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);
    name = sprintf('\\textbf{Mode %d}\n($\\lambda_0$ = %.3g nm, $\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$)', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows,nb_col_mode,u)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, intensities_oseen(:, res_start)+intensities_oseen(:, res_start+1));
    hold on
    h_fl = fill([0 L_L/(L+2*L_L) L_L/(L+2*L_L) 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fl = fill([(L+L_L)/(L+2*L_L) 1 1 (L+L_L)/(L+2*L_L) ], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_fr = fill([L_L/(L+2*L_L) (L+L_L)/(L+2*L_L) (L+L_L)/(L+2*L_L) L_L/(L+2*L_L)], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h1 = plot(X_intensities, intensities_oseen(:, res_start)+intensities_oseen(:, res_start+1));
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows,nb_col_mode,[u+1, nb_col_mode*n_rows]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h_fl;h_fr],'I','LH medium', 'RH medium', 'interpreter', 'latex');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'hybrid_cavity', 'intensity_distribution_oseen', 'eps', 'epsc')

