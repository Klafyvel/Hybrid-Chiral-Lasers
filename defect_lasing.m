% defect_cavity.m
%
% This file is used to simulate an defect chiral cavity with a defect
%
% The cavity looks like this
%
%        interface 1       interface LL           interface 2
% isotropic  | LHCM             | LHCM               | isotropic 
% medium     | n_bar, delta_n   | n_bar_L, delta_n_L | medium
% n_1        | p                | psi, p             | n_2
%             <------ L/2------> <-------- L/2------> 
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
psi = pi/2;                % Rotation of the RH medium
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

% Positions for scanning intensity in LH media
intensity_scan_pos = linspace(0,L/2,400);

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

% In order to improve speed, pre-allocate the array
lasing_coeff = zeros(length(gain_range), length(detuning_range));
% The cavity matrices will be stored here
cavities = zeros(4,4,length(gain_range), length(detuning_range));

% In order to improve speed, pre-allocate the array
lasing_coeff_oseen = zeros(length(gain_range), length(detuning_range));
% The cavity matrices will be stored here
cavities_oseen = zeros(4,4,length(gain_range), length(detuning_range));

for j=1:length(detuning_range)
    for l=1:length(gain_range)
        % Matrix for the first LH media
        M_1 = cwt_from_reduced_variables(p, 0, L/2, k(l,j), kappa, deltak(l,j), lhm);
        % Matrix for the second LH media
        M_2 = cwt_from_reduced_variables(p, psi, L/2, k(l,j), kappa, deltak(l,j), lhm);
        % Matrix for interface 1 (isotrope -> medium 1).
        M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, 0, 0, k(l,j), kappa, n1, k_0(l,j), lhm)^-1;
        % Matrix for interface LL (medium 1 -> medium 2). 
        M_interface_LL = interface_chiral_same_handedness_from_reduced_variables(p, kappa, n_bar(l,j),0, p, kappa, n_bar(l,j), psi, L/2, k_0(l,j),lhm);
        % Matrix for interface 2 (LH -> isotrope)
        M_interface_2 = interface_chiral_to_isotrope_from_reduced_variables(p, psi, L/2, k(l,j), kappa, n2, k_0(l,j), lhm);

        M_cavity = M_interface_2*M_2*M_interface_LL*M_1*M_interface_1;
        lasing_coeff(l,j) = abs(inv(det(M_cavity(3:4,3:4))));
        cavities(:,:,l,j) = M_cavity;

        % Oseen method
        % Transfer matrix inside of the medium
        M_1 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), 0, L/2, k_0(l,j));
        M_2 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), psi, L/2, k_0(l,j));
        % Transfer matrix, accounting for the interfaces, in circular basis.
        M_oseen = T^-1*M_int_2^-1*M_2*M_1*M_int_1*T;
        % Here I use the boundary condition that no field is coming from outside
        lasing_coeff_oseen(l,j) = abs(inv(det(M_oseen(3:4,3:4))));
        cavities_oseen(:,:,l,j) = M_oseen;
    end
end

% %% Same but for the animation.
% % In order to improve speed, pre-allocate the array
% lasing_coeff_anim = zeros(length(psis), length(gain_range), length(detuning_range));
% % The cavity matrices will be stored here
% cavities_anim = zeros(4,4,length(psis), length(gain_range), length(detuning_range));
% 
% for m=1:length(psis)
%     for j=1:length(detuning_range)
%         for l=1:length(gain_range)
%             % Matrix for the first LH media
%             M_1 = cwt_from_reduced_variables(p, 0, L, k(l,j), kappa, deltak(l,j), lhm);
%             % Matrix for the second LH media
%             M_2 = cwt_from_reduced_variables(p, psis(m), L, k(l,j), kappa, deltak(l,j), lhm);
%             % Matrix for interface 1 (isotrope -> medium 1).
%             M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, 0, 0, k(l,j), kappa, n1, k_0(l,j), lhm)^-1;
%             % Matrix for interface LL (medium 1 -> medium 2). 
%             M_interface_LL = interface_chiral_same_handedness_from_reduced_variables(p, kappa, n_bar(l,j),0, p, kappa, n_bar(l,j), psis(m), L, k_0(l,j),lhm);
%             % Matrix for interface 2 (LH -> isotrope)
%             M_interface_2 = interface_chiral_to_isotrope_from_reduced_variables(p, 0, L, k(l,j), kappa, n2, k_0(l,j), lhm);
%     
%             M_cavity = M_interface_2*M_2*M_interface_LL*M_1*M_interface_1;
%             lasing_coeff_anim(m,l,j) = abs(inv(det(M_cavity(3:4,3:4))));
%             cavities_anim(:,:,m,l,j) = M_cavity;
%         end
%     end
% end

%% Finding lasing loci by looking for peaks in the log of the lasing coefficient
% The modes are identified by peaks on the logarithm of the invert of the 
% determinant. 
coef_log = log(lasing_coeff);
% row and col will store the rows and columns where lasing modes are found
[row, col] = find_max(coef_log, peak_radius/step_detuning, peak_radius/step_gain);

coef_log_oseen = log(lasing_coeff_oseen);
% row_oseen and col_oseen will store the rows and columns where lasing modes are 
% found with oseen transformation 
[row_oseen, col_oseen] = find_max(coef_log_oseen, peak_radius/step_detuning, peak_radius/step_gain);

% coef_log_anim = log(lasing_coeff_anim);
% for m=1:length(psis)
%     [row_anim_m, col_anim_m] = find_max(squeeze(coef_log_anim(m,:,:)), peak_radius/step_detuning, peak_radius/step_gain);
%     row_anim{m} = row_anim_m;
%     col_anim{m} = col_anim_m;
% end


%% Calculate the field at z=0 for each mode 

% R, theta and epsilon define the ellipse, lcp the handedness (true if left-handed) and modes are the eigen vectors, i.e. the output modes.
[R, theta, epsilon, lcp, modes] = find_modes(row,col,cavities);

% Same but for the Oseen method.
[R_oseen, theta_oseen, epsilon_oseen, lcp_oseen, modes_oseen] = find_modes(row_oseen,col_oseen,cavities_oseen);

% for m=1:length(psis)
%     [R_m, theta_m, epsilon_m, lcp_m, modes_m] = find_modes(row_anim{m},col_anim{m},squeeze(cavities_anim(:,:,m,:,:)));
%     R_anim{m} = R_m;
%     theta_anim{m} = theta_m;
%     epsilon_anim{m} = epsilon_m;
%     lcp_anim{m} = lcp_m;
%     modes_anim{m} = modes_m;
% end

%% Calculate the intensity distribution

% Number of scan positions in the cavity
nb_scan_pos = length(intensity_scan_pos);

% the intensities will be stored here
intensities = zeros(2*nb_scan_pos, 4*length(row));

for u=1:length(row)
    l = row(u);
    j = col(u);

    m = modes(:,u);

    % Field at z=0⁻
    F_0 = [0;0;m(1);m(2)];

    % Transfer matrix at interface 1
    M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, 0, 0, k(l,j), kappa, n1, k_0(l,j), lhm)^-1;
    % Matrix for interface LL (medium 1 -> medium 2). 
    M_interface_LL = interface_chiral_same_handedness_from_reduced_variables(p, kappa, n_bar(l,j),0, p, kappa, n_bar(l,j), psi, L/2, k_0(l,j),lhm);

    % Shortcuts to store the results of intensity scan at the right place
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan the first medium
    for v=1:nb_scan_pos
        % length of the partial slab
        len = intensity_scan_pos(v);

        % Corresponding matrix
        M_1 = cwt_from_reduced_variables(p, 0, len, k(l,j), kappa, deltak(l,j), lhm);
        % Field at z
        F = M_1 * M_interface_1 * F_0;

        % Intensity at z
        intensities(v,res_start:res_stop) = abs(F).^2;
    end
    % Shortcut to store the matrix corresponding to the first medium
    M = M_interface_LL * cwt_from_reduced_variables(p, 0, L/2, k(l,j), kappa, deltak(l,j), lhm) * M_interface_1;

    % Scan the second medium
    for v=1:nb_scan_pos
        % length of the partial slab
        len = intensity_scan_pos(v);

        % Corresponding matrix
        M_2 = cwt_from_reduced_variables(p, psi, len, k(l,j), kappa, deltak(l,j), lhm);
        % Field at z
        F = M_2*(M * F_0);
        % Intensity at z
        intensities(v+nb_scan_pos,res_start:res_stop) = abs(F).^2;
    end
end

% %% Calculate the intensity distribution (animation)
% 
% % the intensities will be stored here
% intensities_anim = cell(length(psis));
% 
% for m=1:length(psis)
%     intensities_anim{m} = zeros(2*nb_scan_pos, 4*length(row_anim{m}));
%     for u=1:length(row_anim{m})
%         l = row_anim{m}(u);
%         j = col_anim{m}(u);
%     
%         mode = modes_anim{m}(:,u);
%     
%         % Field at z=0⁻
%         F_0 = [0;0;mode(1);mode(2)];
%     
%         % Transfer matrix at interface 1
%         M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, 0, 0, k(l,j), kappa, n1, k_0(l,j), lhm)^-1;
%         % Matrix for interface LL (medium 1 -> medium 2). 
%         M_interface_LL = interface_chiral_same_handedness_from_reduced_variables(p, kappa, n_bar(l,j),0, p, kappa, n_bar(l,j), psis(m), L, k_0(l,j),lhm);
%     
%         % Shortcuts to store the results of intensity scan at the right place
%         res_start = 4*(u-1)+1;
%         res_stop = 4*u;
%     
%         % Scan the first medium
%         for v=1:nb_scan_pos
%             % length of the partial slab
%             len = intensity_scan_pos(v);
%     
%             % Corresponding matrix
%             M_1 = cwt_from_reduced_variables(p, 0, len, k(l,j), kappa, deltak(l,j), lhm);
%             % Field at z
%             F = M_1 * M_interface_1 * F_0;
%     
%             % Intensity at z
%             intensities_anim{m}(v,res_start:res_stop) = abs(F).^2;
%         end
%         % Shortcut to store the matrix corresponding to the first medium
%         M = M_interface_LL * cwt_from_reduced_variables(p, 0, L, k(l,j), kappa, deltak(l,j), lhm) * M_interface_1;
%     
%         % Scan the second medium
%         for v=1:nb_scan_pos
%             % length of the partial slab
%             len = intensity_scan_pos(v);
%     
%             % Corresponding matrix
%             M_2 = cwt_from_reduced_variables(p, psis(m), len, k(l,j), kappa, deltak(l,j), lhm);
%             % Field at z
%             F = M_2*(M * F_0);
%             % Intensity at z
%             intensities_anim{m}(v+nb_scan_pos,res_start:res_stop) = abs(F).^2;
%         end
%     end
% end
%% Scan through the cavity with the Oseen method to check SVEA

% Pre-allocate the array for storing intensities
intensities_oseen = zeros(2*nb_scan_pos, 4*length(row_oseen));

for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);

    m = modes_oseen(:,u);

    % Field at z=0⁻ (electromagnetic basis)
    F_0 = M_int_1*T*[0;0;m(1);m(2)];

    % Shorthand to access the array of intensities at the right position
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan through first medium
    for v=1:nb_scan_pos
        % position in the medium, -> length of the virtual slab 
        len = intensity_scan_pos(v);
        % Corresponding matrix
        M_1 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), 0, len, k_0(l,j));
        % Intensity at z
        intensities_oseen(v,res_start:res_stop) = abs(M_1*F_0).^2;
    end
    % Shorthand for the first medium
    M = oseen(p, epsilon_a(l,j), epsilon_b(l,j), 0, L/2, k_0(l,j));

    % Scan through second medium
    for v=1:nb_scan_pos
        % Position in RH medium
        len = intensity_scan_pos(v);
        % Corresponding matrix
        M_2 = oseen(p, epsilon_a(l,j), epsilon_b(l,j), psi, len, k_0(l,j));
        % Intensity at z
        intensities_oseen(v+nb_scan_pos,res_start:res_stop) = abs(M_2*M*F_0).^2;
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

savefig(gcf, 'defect_cavity', 'lasing', 'png')

% Plot for the report
figure('visible', 'off')
surf(X,Y,coef_log, 'EdgeColor', 'none')
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
zlabel('$-\log(\mathrm{det})$', 'interpreter', 'latex')
savefig(gcf, 'defect_cavity', 'surface', 'eps', 'epsc')
figure('visible', 'off')
contour(X,Y,coef_log, 30)
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
for j=1:length(row)
    x = X(row(j), col(j));
    y = Y(row(j), col(j));
    s = sprintf('$\\leftarrow %.3g %+.3g i$', x, y);
    text(x, y, s, 'interpreter', 'latex')
end
savefig(gcf, 'defect_cavity', 'contour', 'eps', 'epsc')


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

savefig(gcf, 'defect_cavity', 'lasing_oseen', 'png')

% Plot for the report
figure('visible', 'off')
surf(X,Y,coef_log_oseen, 'EdgeColor', 'none')
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
zlabel('$-\log(\mathrm{det})$', 'interpreter', 'latex')
savefig(gcf, 'defect_cavity', 'surface_oseen', 'eps', 'epsc')
figure('visible', 'off')
contour(X,Y,coef_log_oseen, 30)
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
for j=1:length(row_oseen)
    x = X(row_oseen(j), col_oseen(j));
    y = Y(row_oseen(j), col_oseen(j));
    s = sprintf('$\\leftarrow %.3g %+.3g i$', x, y);
    text(x, y, s, 'interpreter', 'latex')
end
savefig(gcf, 'defect_cavity', 'contour_oseen', 'eps', 'epsc')

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

X_intensities = cat(2, intensity_scan_pos, intensity_scan_pos+L/2) ./ L;

for u=1:length(row)
    l = row(u);
    j = col(u);
    name = sprintf('Mode %d\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows+2,nb_col_mode,u+2*nb_col_mode)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, res_start), '-','LineWidth',2);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    h_f1 = fill([0 .5 .5 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_f2 = fill([.5 1 1 .5], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
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
lgd = legend(hL,[h1;h2;h3;h4;h5;h_f1;h_f2],'I','I_L^+','I_R^+','I_L^-','I_R^-','Medium 1', 'Medium 2');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'defect_cavity', 'intensity_distribution', 'png')

% plot for the report
figure('visible', 'off')
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)))
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
daspect([1 1 1])
savefig(gcf, 'defect_cavity', 'modes_found', 'eps', 'epsc')

figure('visible', 'off', 'Position', [0 0 width n_rows*250])
for u=1:length(row)
    l = row(u);
    j = col(u);
    name = sprintf('Mode %d\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows,nb_col_mode,u)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities(:,res_start:res_stop), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, res_start), '-','LineWidth',2);
    h3 = plot(X_intensities, intensities(:, (res_start+1)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (res_start+2)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:, (res_start+3)), '-','LineWidth',1);
    h_f1 = fill([0 .5 .5 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_f2 = fill([.5 1 1 .5], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
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
hL = subplot(n_rows,nb_col_mode,[u+1 nb_col_mode*n_rows]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h2;h3;h4;h5;h_f1;h_f2],'$I$','$I_L^+$','$I_R^+$','$I_L^-$','$I_R^-$','Medium 1', 'Medium 2','NumColumns',3, 'interpreter', 'latex');
set(lgd,'position',poshL);
axis(hL,'off');
savefig(gcf, 'defect_cavity', 'intensity_distribution', 'eps', 'epsc')

%% Plot field modulus distribution from the Oseen method

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

for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);
    name = sprintf('Mode %d\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows+2,nb_col_mode,u+2*nb_col_mode)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities_oseen(:, res_start:res_stop),2));
    hold on
    h_f1 = fill([0 .5 .5 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_f2 = fill([.5 1 1 .5], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h1 = plot(X_intensities, sum(intensities_oseen(:, res_start:res_stop),2));
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows+2,nb_col_mode,[u+1+2*nb_col_mode nb_col_mode*(n_rows+2)]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h_f1;h_f2],'I','Medium 1', 'Medium 2');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'defect_cavity', 'intensity_distribution_oseen', 'png')

% plot for the report

figure('visible', 'off')
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
xlabel('Detuning', 'interpreter', 'latex')
ylabel('Gain', 'interpreter', 'latex')
daspect([1 1 1])
savefig(gcf, 'defect_cavity', 'modes_found_oseen', 'eps', 'epsc')

figure('visible', 'off', 'Position', [0 0 width n_rows*250])
for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);
    name = sprintf('Mode %d\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows,nb_col_mode,u)
    res_start = 4*(u-1)+1;
    res_stop = 4*u;
    h1 = plot(X_intensities, sum(intensities_oseen(:, res_start:res_stop),2));
    hold on
    h_f1 = fill([0 .5 .5 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
    h_f2 = fill([.5 1 1 .5], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
    h1 = plot(X_intensities, sum(intensities_oseen(:, res_start:res_stop),2));
    title(name, 'interpreter', 'latex')
    xlabel('z/length', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows,nb_col_mode,[u+1, nb_col_mode*n_rows]);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h_f1;h_f2],'I','Medium 1', 'Medium 2', 'interpreter', 'latex');
set(lgd,'position',poshL);
axis(hL,'off');

savefig(gcf, 'defect_cavity', 'intensity_distribution_oseen', 'eps', 'epsc')
% %% NOW IT'S SHOW TIME
% close all;
% figure
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% loops = length(psis);
% % initialise the frames
% frames(loops) = getframe(gcf);
% 
% v = VideoWriter('plots/defect.avi');
% open(v);
% 
% for m=1:loops
%     n_rows = ceil(length(row_anim{m}) / nb_col_mode);
%     if mod(length(row_anim{m}), nb_col_mode) == 0
%         n_rows = n_rows + 1;
%     end
%     figure('Name', 'Intensity distribution for found modes', 'Position', [0 0 1200 (n_rows+2)*250], 'visible', 'off')
%     subplot(n_rows+2,nb_col_mode,[1 2*nb_col_mode])
%     plot_ellipse(epsilon_anim{m},theta_anim{m},lcp_anim{m},X(row_anim{m}(1),col_anim{m}),Y(row_anim{m},col_anim{m}(1)))
%     xlim([-1 7])
%     ylim([-1 3])
%     name = sprintf('Output modes for cavity with a defect\n$\\psi=%.3g^{\\:\\circ}$', psis(m)*180/pi);
%     title(name, 'interpreter', 'latex')
%     xlabel('Detuning', 'interpreter', 'latex')
%     ylabel('Gain', 'interpreter', 'latex')
%     daspect([1 1 1])
%     
%     X_intensities = cat(2, intensity_scan_pos, intensity_scan_pos+L) ./ (2*L);
%     
%     for u=1:length(row_anim{m})
%         l = row_anim{m}(u);
%         j = col_anim{m}(u);
%         name = sprintf('Mode %d\n($\\lambda_0$ = %.3g nm, $\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$)', u, 2*pi/k_0(l,j)*1e9, theta_anim{m}(u)*180/pi, epsilon_anim{m}(u)*180/pi);
%         subplot(n_rows+2,nb_col_mode,u+2*nb_col_mode)
%         res_start = 4*(u-1)+1;
%         res_stop = 4*u;
%         h1 = plot(X_intensities, sum(intensities_anim{m}(:,res_start:res_stop), 2),'LineWidth',2);
%         hold on
%         h2 = plot(X_intensities, intensities_anim{m}(:, res_start), '-','LineWidth',2);
%         h3 = plot(X_intensities, intensities_anim{m}(:, (res_start+1)), '-','LineWidth',1);
%         h4 = plot(X_intensities, intensities_anim{m}(:, (res_start+2)), '-','LineWidth',1);
%         h5 = plot(X_intensities, intensities_anim{m}(:, (res_start+3)), '-','LineWidth',1);
%         h_f1 = fill([0 .5 .5 0], [max(ylim) max(ylim) min(ylim) min(ylim)],'m','FaceAlpha',0.1, 'LineStyle', 'none');
%         h_f2 = fill([.5 1 1 .5], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','FaceAlpha',0.1, 'LineStyle', 'none');
%         h1 = plot(X_intensities, sum(intensities_anim{m}(:,res_start:res_stop), 2),'-','LineWidth',2,'Color', [247 216 79]/255);
%         h2 = plot(X_intensities, intensities_anim{m}(:, res_start), 'r-','LineWidth',1);
%         h3 = plot(X_intensities, intensities_anim{m}(:, (res_start+1)), 'b-','LineWidth',1);
%         h4 = plot(X_intensities, intensities_anim{m}(:, (res_start+2)), '-','LineWidth',1);
%         h5 = plot(X_intensities, intensities_anim{m}(:, (res_start+3)), '-','LineWidth',1);
%         title(name, 'interpreter', 'latex')
%         xlabel('z/length', 'interpreter', 'latex')
%         ylabel('I (a.u.)', 'interpreter', 'latex')
%         hold off
%     end
%     hL = subplot(n_rows+2,nb_col_mode, u+2*nb_col_mode+1);
%     poshL = get(hL,'position');
%     lgd = legend(hL,[h1;h2;h3;h4;h5;h_f1;h_f2],'I','I_L^+','I_R^+','I_L^-','I_R^-','Medium 1', 'Medium 2');
%     set(lgd,'position',poshL);
%     axis(hL,'off');
% 
%     frame = getframe(gcf);
%     frames(m) = frame;
%     writeVideo(v,frame);
% end
% 
% close(v);
% % play the movie @ 24 fps
% movie(frames,2, 24);
% 
% %% Same but with lasing_coeff
% 
% close all;
% clear frames;
% figure
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% loops = length(psis);
% % initialise the frames
% frames(loops) = getframe(gcf);
% 
% v = VideoWriter('plots/lasing_coeff_defect.avi');
% open(v);
% for m=1:length(psis)
%     figure('Position', [0 0 1200 600], 'visible', 'off')
%     subplot(121)
%     surf(X,Y,squeeze(coef_log_anim(m,:,:)), 'EdgeColor', 'none')
%     xlabel('Detuning')
%     ylabel('Gain')
%     zlim([-10 5])
%     ylim([-1 2])
%     name = sprintf('Lasing coefficient for cavity with a defect\n$\\psi=%.3g^{\\:\\circ}$', psis(m)*180/pi);
%     title(name, 'interpreter', 'latex')
%     subplot(122)
%     contour(X,Y,squeeze(coef_log_anim(m,:,:)), 30)
%     xlabel('Detuning')
%     ylabel('Gain')
%     ylim([-1 2])
%     for j=1:length(row_anim{m})
%         x = X(row_anim{m}(j), col_anim{m}(j));
%         y = Y(row_anim{m}(j), col_anim{m}(j));
%         s = sprintf('< %.3g %+.3g i', x, y);
%         text(x, y, s)
%     end
%     frame = getframe(gcf);
%     frames(m) = frame;
%     writeVideo(v,frame);
% end
% close(v);
% % play the movie @ 24 fps
% movie(frames,2, 24);
% 
% 
