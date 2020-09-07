% simple_cavity.m
%
% This file simulates a simple chiral cavity.
% 
% Plot the invert of the determinant of the matrix giving the transformation for a round trib in a chiral cavity.
close all;
clear all;
%% Data of the simulation
Lp = 300e-9;            % "wavelength" of the chiral medium in m

L = 20*Lp;              % Length of the cavity in m
kappaL = 4;             % kappa * L where kappa is the coupling constant
psi = 0;                % Rotation in x-y plane of the medium

re_n_bar = 1.7690;      % Real part of average refractive index of the medium
n1 = 1;                 % Refractive index of medium 1
n2 = n1;                % Refractive index of medium 2

max_detuning = 12;      % maximum detuning
max_gain = 2.5,         % maximum gain
step_gain = 0.02;       % Step for the range of gains explored
step_detuning = 0.02;   % Step for the range of detunings explored

scan_pos = linspace(0,L,300); % position scanned for each modes to determine intensity

peak_radius = 0.2;      % Expected radius of lasing peaks in the lasing plot. 
                        % Reduce if you think some peaks are missing, increase 
                        % if you think there are false positives.

%% Useful intermediary variables
p = 2*pi/Lp;            % Periodicity of the chiral medium
kappa = kappaL / L;     % Coupling constant
I = eye(4);             % Identity matrix

gain_range = -max_gain*0.1:step_gain:max_gain; % Gain range explored (-Im(\delta k \times L / 2))
detuning_range = 0:step_detuning:max_detuning; % Detuning range explored (Re(\delta k \times L / 2))

r_1 = (re_n_bar - n1) / (re_n_bar + n1); % reflection coefficient on interface 1 for Toph and McCall paper
r_2 = (re_n_bar - n2) / (re_n_bar + n2); % reflection coefficient on interface 1 for Toph and McCall paper
t_1 = 2*re_n_bar / (re_n_bar + n1); % transmission coefficient on interface 1 for Toph and McCall paper
t_2 = 2*re_n_bar / (re_n_bar + n2); % transmission coefficient on interface 1 for Toph and McCall paper

[X,Y] = meshgrid(detuning_range,gain_range); % X and Y mesh for detning and gain ranges
deltak = 2.*X ./ L - 2.*i.*Y./L; % detuning
k = deltak ./ 2 + p;    % Wavevector in medium
k_0 = real(k) ./ re_n_bar;       % Wavevector in vacuum
n_bar = k ./ k_0;

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
%% Compute invert of determinant of the reflection matrices
% The determinant of the matrix giving a round trip in the cavity will be stored here.
lasing_coeff = zeros(length(gain_range), length(detuning_range));
lasing_coeff_topf = zeros(length(gain_range), length(detuning_range));
lasing_coeff_oseen = zeros(length(gain_range), length(detuning_range));
% lasing_coeff_dummy = zeros(length(gain_range), length(detuning_range));
% The cavity matrices will be stored here
cavities = zeros(4,4,length(gain_range), length(detuning_range));
cavities_oseen = zeros(4,4,length(gain_range), length(detuning_range));
% cavities_dummy = zeros(4,4,length(gain_range), length(detuning_range));

for j=1:length(detuning_range)
    for l=1:length(gain_range)
        % Transfer matrix in the medium
        M_medium = cwt_from_reduced_variables(p, psi, L, k(l,j), kappa, deltak(l,j));
        % Transfer matrix for interfaces 1 and 2
        M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, psi, 0, k(l,j), kappa, n1, k_0(l,j));
        M_interface_2 = interface_chiral_to_isotrope_from_reduced_variables(p, psi, L, k(l,j), kappa, n2, k_0(l,j));
        cavities(:,:,l,j) = M_interface_2*M_medium*M_interface_1^-1;

        % Below there is a calculation of the lasing coefficient with what is called 
        % 'The first method' in the report. But for convenience I decided to use the same 
        % method everywhere, a.k.a. 'second method'.
        % Reflectivity matrix for the medium, i.e. 
        % [E_L^+(L^-), E_R^+(L^-), E_L^-(0^+), E_R^-(0^+)]^T = r[E_L^+(0^+), E_R^+(0^+), E_L^-(L^-), E_R^-(L^-)]^T
        % reflectivity_medium = partial_inverse(M_medium)^-1;

        % % Reflectivity matrices for both interfaces (2x2 matrices)
        % R_1 = reflection_cwt(M_interface_1);
        % R_2 = reflection_cwt(M_interface_2);

        % % Boundary condition is given by the reflections at the interfaces
        % boundary_condition = [
        %     zeros(2) R_1;
        %     R_2 zeros(2);
        % ];

        % % Matrix M gives the field just before the interfaces after one round trip.
        % M = boundary_condition*reflectivity_medium;
        % % In order to achieve lasing, M should be identity
        % lasing_coeff(l,j) = inv(abs(det(M - I)));

        
        lasing_coeff(l,j) = abs(inv(det(cavities(3:4,3:4,l,j))));


        % Topf and McCall paper
        Delta = sqrt((kappa)^2-(deltak(l,j)/2)^2);
        P_plus = exp(i*p*L) * (cosh(Delta*L)+i*deltak(l,j)/2/Delta*sinh(Delta*L));
        P_minus = exp(-i*p*L) * (cosh(Delta*L)-i*deltak(l,j)/2/Delta*sinh(Delta*L));
        lasing_coeff_topf(l,j) = abs(inv((r_1*r_2*exp(i*k(l,j)*L))^2*P_plus-2*r_1*r_2*exp(i*k(l,j)*L)+P_minus));

        % Oseen method
        % Transfer matrix inside of the medium
        M_oseen_medium = oseen(p, epsilon_a(l,j), epsilon_b(l,j), 0, L, k_0(l,j));
        % Transfer matrix, accounting for the interfaces, in circular basis.
        M_oseen = T^-1*M_int_2^-1*M_oseen_medium*M_int_1*T;
        cavities_oseen(:,:,l,j) = M_oseen;
        % Here I use the boundary condition that no field is coming from outside
        lasing_coeff_oseen(l,j) = abs(inv(det(M_oseen(3:4,3:4))));

        % Below there is a 'dummy' cwt, meant to check the results given with the interface matrices
        % For CWT is correct.
        % % Dummy cwt
        % M_L = cwt_circular_to_emag_basis_from_reduced_variables(p,psi,kappa,n_bar(l,j),0,k_0(l,j));
        % dummy(:,:,l,j) = inv(M_L)*M_int_1*T;
        % M_R = cwt_circular_to_emag_basis_from_reduced_variables(p,psi,kappa,n_bar(l,j),L,k_0(l,j));
        % M_dummy = inv(T)*inv(M_int_2)*M_R * M_medium * inv(M_L)*M_int_1*T;
        % cavities_dummy(:,:,l,j) = M_dummy;
        % lasing_coeff_dummy(l,j) = abs(inv(det(M_dummy(3:4,3:4))));
    end
end

%% Finding lasing loci using the plots

% The modes are identified by pikes on the logarithm of the invert of the 
% determinant. 
coef_log = log(lasing_coeff);
% row and col will store the rows and columns where lasing modes are found
[row, col] = find_max(coef_log, peak_radius/step_detuning, peak_radius/step_gain);

coef_log_topf = log(lasing_coeff_topf);
% row_topf and col_topf will store the rows and columns where lasing modes are found
% with Topf and McCall approach
[row_topf, col_topf] = find_max(coef_log_topf, peak_radius/step_detuning, peak_radius/step_gain);

coef_log_oseen = log(lasing_coeff_oseen);
% row_oseen and col_oseen will store the rows and columns where lasing modes are 
% found with oseen transformation 
[row_oseen, col_oseen] = find_max(coef_log_oseen, peak_radius/step_detuning, peak_radius/step_gain);

% coef_log_dummy = log(lasing_coeff_dummy);
% % row_dummy and col_dummy will store the rows and columns where lasing modes are 
% % found with oseen transformation 
% [row_dummy, col_dummy] = find_max(coef_log_dummy, peak_radius/step_detuning, peak_radius/step_gain);

%% Calculate the field at z=0 for each mode

% R, theta and epsilon define the ellipse, lcp the handedness (true if left-handed) and modes are the eigen vectors, i.e. the output modes.
[R, theta, epsilon, lcp, modes] = find_modes(row,col,cavities);

% Same but for the Oseen method.
[R_oseen, theta_oseen, epsilon_oseen, lcp_oseen, modes_oseen] = find_modes(row_oseen,col_oseen,cavities_oseen);


%% Calculate intensity distribution

% the intensities will be stored here
intensities = zeros(length(scan_pos), 4*length(row));

for u=1:length(row)
    l = row(u);
    j = col(u);

    m = modes(:,u);
    % Field at z=0⁻
    F_0 = [0;0;m(1);m(2)];

    % Transfer matrix for interfaces 1 and 2
    M_interface_1 = interface_chiral_to_isotrope_from_reduced_variables(p, psi, 0, k(l,j), kappa, n1, k_0(l,j))^-1;

    % Shortcuts to store the results of intensity scan at the right place
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan through the medium
    for v=1:length(scan_pos)
        pos = scan_pos(v);
        M_partial = cwt_from_reduced_variables(p, psi, pos, k(l,j), kappa, deltak(l,j));
        field_at_pos = M_partial * M_interface_1 * F_0;
        intensities(v,res_start:res_stop) = abs(field_at_pos) .^ 2;
    end
end

%% Calculate intensity distribution (Oseen)

% the intensities will be stored here
intensities_oseen = zeros(length(scan_pos), 4*length(row));

for u=1:length(row_oseen)
    l = row_oseen(u);
    j = col_oseen(u);

    m = modes_oseen(:,u);
    % Field at z=0⁻
    F_0 = M_int_1*T*[0;0;m(1);m(2)];

    % Shortcuts to store the results of intensity scan at the right place
    res_start = 4*(u-1)+1;
    res_stop = 4*u;

    % Scan through the medium
    for v=1:length(scan_pos)
        pos = scan_pos(v);
        M_partial = oseen(p, epsilon_a(l,j), epsilon_b(l,j), psi, pos, k_0(l,j));
        field_at_pos = M_partial * F_0;
        intensities_oseen(v,res_start:res_stop) = abs(field_at_pos) .^ 2;
    end
end

%% %%%%%%%%% BEYOND THAT POINT, THERE ARE ONLY PLOTS. %%%%%%%%% %%

%% Plot the surface and contour


figure('Position', [0 0 1200 1600])
subplot(321)
surf(X,Y,coef_log, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Developped CWT')
subplot(322)
contour(X,Y,coef_log, 30)
xlabel('Detuning')
ylabel('Gain')
title('Developped CWT')
for j=1:length(row)
    s = sprintf('< %.3g %+.3g i', X(row(j), col(j)), Y(row(j), col(j)));
    text(X(row(j), col(j)),Y(row(j), col(j)),s)
end
subplot(323)
surf(X,Y,coef_log_topf, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Topf and McCall 2014 paper')
subplot(324)
contour(X,Y,coef_log_topf, 30)
xlabel('Detuning')
ylabel('Gain')
title('Topf and McCall 2014 paper')
for j=1:length(row_topf)
    s = sprintf('< %.3g %+.3g i', X(row_topf(j), col_topf(j)),Y(row_topf(j), col_topf(j)));
    text(X(row_topf(j), col_topf(j)),Y(row_topf(j), col_topf(j)),s)
end
subplot(325)
surf(X,Y,coef_log_oseen, 'EdgeColor', 'none')
xlabel('Detuning')
ylabel('Gain')
title('Oseen transformation')
subplot(326)
contour(X,Y,coef_log_oseen, 30)
xlabel('Detuning')
ylabel('Gain')
title('Oseen transformation')
for j=1:length(row_oseen)
    s = sprintf('< %.3g %+.3g i', X(row_oseen(j), col_oseen(j)),Y(row_oseen(j), col_oseen(j)));
    text(X(row_oseen(j), col_oseen(j)),Y(row_oseen(j), col_oseen(j)),s)
end
% subplot(427)
% surf(X,Y,coef_log_dummy, 'EdgeColor', 'none')
% xlabel('Detuning')
% ylabel('Gain')
% title('Dummy transformation')
% subplot(428)
% contour(X,Y,coef_log_dummy, 30)
% xlabel('Detuning')
% ylabel('Gain')
% title('Dummy transformation')
% for j=1:length(row_dummy)
%     s = sprintf('< %.3g %+.3g i', X(row_dummy(j), col_dummy(j)),Y(row_dummy(j), col_dummy(j)));
%     text(X(row_dummy(j), col_dummy(j)),Y(row_dummy(j), col_dummy(j)),s)
% end


savefig(gcf, 'simple_cavity', 'lasing', 'png')

%% For the report, it is more convenient to have separated figures
figure('visible','off')
surf(X,Y,coef_log, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'simple_cavity', 'surface', 'eps', 'epsc')
figure('visible','off')
contour(X,Y,coef_log, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row)
    s = sprintf('$\\leftarrow %.3g %+.3g i$', X(row(j), col(j)), Y(row(j), col(j)));
    text(X(row(j), col(j)),Y(row(j), col(j)),s,'interpreter','latex')
end
savefig(gcf, 'simple_cavity', 'contour', 'eps', 'epsc')
figure('visible','off')
surf(X,Y,coef_log_topf, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'simple_cavity', 'surface_topf', 'eps', 'epsc')
figure('visible','off')
contour(X,Y,coef_log_topf, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row_topf)
    s = sprintf('$\\leftarrow %.3g %+.3g i$', X(row_topf(j), col_topf(j)),Y(row_topf(j), col_topf(j)));
    text(X(row_topf(j), col_topf(j)),Y(row_topf(j), col_topf(j)),s,'interpreter','latex')
end
savefig(gcf, 'simple_cavity', 'contour_topf', 'eps', 'epsc')

figure('visible','off')
surf(X,Y,coef_log_oseen, 'EdgeColor', 'none')
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
zlabel('$-\log(\mathrm{det})$','interpreter','latex')
savefig(gcf, 'simple_cavity', 'surface_oseen', 'eps', 'epsc')
figure('visible','off')
contour(X,Y,coef_log_oseen, 30)
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
for j=1:length(row_oseen)
    s = sprintf('$\\leftarrow %.3g %+.3g i$', X(row_oseen(j), col_oseen(j)),Y(row_oseen(j), col_oseen(j)));
    text(X(row_oseen(j), col_oseen(j)),Y(row_oseen(j), col_oseen(j)),s,'interpreter','latex')
end
savefig(gcf, 'simple_cavity', 'contour_oseen', 'eps', 'epsc')

%% Plot modes intensity distribution

n_rows = round(length(row) / 4);
if mod(length(row), 4) == 0
    n_rows = n_rows + 1;
end
figure('Name', 'Analysis of found modes','Position', [0 0 1200 1000])
subplot(n_rows+2,4,[1 8])
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)))
title('Identification of output modes, developped CWT')
xlabel('Detuning')
ylabel('Gain')
daspect([1 1 1])

X_intensities = linspace(0,1,length(scan_pos));
for u=1:length(row)
    name = sprintf('\\textbf{Mode %d}\n$\\lambda_0$ = %.3g nm\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows+2,4,u+8)
    h1 = plot(X_intensities, sum(intensities(:,(4*(u-1)+1):4*u), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, (4*(u-1)+1)), '-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (4*(u-1)+2)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (4*(u-1)+3)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:,(4*(u-1)+4)), '-','LineWidth',1);
    title(name, 'interpreter', 'latex')
    xlabel('z/L', 'interpreter', 'latex')
    ylabel('Intensity (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows+2,4,u+9);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h2;h3;h4;h5],'I','I_L^+','I_R^+','I_L^-','I_R^-');
set(lgd,'position',poshL);
axis(hL,'off');
savefig(gcf, 'simple_cavity', 'modes_found', 'png')

% For the report, it is more convenient to have separated figures
figure('Name', 'Analysis of found modes', 'visible','off')
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)))
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
daspect([1 1 1])
savefig(gcf, 'simple_cavity', 'modes_found', 'eps', 'epsc')

figure('Name', 'Intensity distribution', 'Position', [0 0 1000 500], 'visible','off')
for u=1:length(row)
    name = sprintf('\\textbf{Mode %d}\n$\\lambda_0$ = %.3g nm\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta(u)*180/pi, epsilon(u)*180/pi);
    subplot(n_rows,4,u)
    h1 = plot(X_intensities, sum(intensities(:,(4*(u-1)+1):4*u), 2),'LineWidth',2);
    hold on
    h2 = plot(X_intensities, intensities(:, (4*(u-1)+1)), '-','LineWidth',1);
    h3 = plot(X_intensities, intensities(:, (4*(u-1)+2)), '-','LineWidth',1);
    h4 = plot(X_intensities, intensities(:, (4*(u-1)+3)), '-','LineWidth',1);
    h5 = plot(X_intensities, intensities(:,(4*(u-1)+4)), '-','LineWidth',1);
    title(name,'interpreter','latex')
    xlabel('z/L','interpreter','latex')
    ylabel('I (a.u.)','interpreter','latex')
    hold off
end
hL = subplot(n_rows,4,u+1);
poshL = get(hL,'position');
lgd = legend(hL,[h1;h2;h3;h4;h5],'$I$','$I_L^+$','$I_R^+$','$I_L^-$','$I_R^-$','interpreter','latex');
set(lgd,'position',poshL);
axis(hL,'off');
savefig(gcf, 'simple_cavity', 'intensity_distribution', 'eps', 'epsc')

%% Plot modes intensity distribution (oseen)

n_rows = round(length(row_oseen) / 4);
if mod(length(row_oseen), 4) == 0
    n_rows = n_rows + 1;
end
figure('Name', 'Analysis of found modes','Position', [0 0 1200 1000])
subplot(n_rows+2,4,[1 8])
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
title('Identification of output modes, Oseen')
xlabel('Detuning')
ylabel('Gain')
daspect([1 1 1])

X_intensities = linspace(0,1,length(scan_pos));
for u=1:length(row_oseen)
    name = sprintf('\\textbf{Mode %d}\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows+2,4,u+8)
    h1 = plot(X_intensities, sum(intensities_oseen(:,(4*(u-1)+1):4*u), 2),'LineWidth',1);
    hold on
    title(name, 'interpreter', 'latex')
    xlabel('z/L', 'interpreter', 'latex')
    ylabel('I (a.u.)', 'interpreter', 'latex')
    hold off
end
hL = subplot(n_rows+2,4,u+9);
poshL = get(hL,'position');
lgd = legend(hL,[h1],'I');
set(lgd,'position',poshL);
axis(hL,'off');
savefig(gcf, 'simple_cavity', 'modes_found_oseen', 'png')

% For the report, it is more convenient to have separated figures
figure('Name', 'Analysis of found modes', 'visible','off')
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
daspect([1 1 1])
savefig(gcf, 'simple_cavity', 'modes_found_oseen', 'eps', 'epsc')

figure('Name', 'Intensity distribution', 'Position', [0 0 1000 500], 'visible','off')
for u=1:length(row_oseen)
    name = sprintf('\\textbf{Mode %d}\n$\\lambda_0$ = %.3g nm,\n$\\theta=%.3g^{\\:\\circ}$, $\\epsilon=%.3g^{\\:\\circ}$', u, 2*pi/k_0(l,j)*1e9, theta_oseen(u)*180/pi, epsilon_oseen(u)*180/pi);
    subplot(n_rows,4,u)
    h1 = plot(X_intensities, sum(intensities_oseen(:,(4*(u-1)+1):4*u), 2),'LineWidth',1);
    hold on
    title(name,'interpreter','latex')
    xlabel('z/L','interpreter','latex')
    ylabel('I (a.u.)','interpreter','latex')
    hold off
end
hL = subplot(n_rows,4,u+1);
poshL = get(hL,'position');
lgd = legend(hL,[h1],'$I$','interpreter','latex');
set(lgd,'position',poshL);
axis(hL,'off');
savefig(gcf, 'simple_cavity', 'intensity_distribution_oseen', 'eps', 'epsc')

% 
% figure('visible','off')
% plot_ellipse(epsilon_topf,theta_topf,lcp_topf,X(row_topf(1),col_topf),Y(row_topf,col_topf(1)))
% xlabel('Detuning','interpreter','latex')
% ylabel('Gain','interpreter','latex')
% daspect([1 1 1])
% savefig(gcf, 'simple_cavity', 'modes_found_topf', 'eps', 'epsc')

%% Plot modes ellipses

% Plot for the report
figure('visible','off')
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)))
xlabel('Detuning','interpreter','latex')
ylabel('Gain','interpreter','latex')
daspect([1 1 1])
savefig(gcf, 'simple_cavity', 'modes_found_oseen', 'eps', 'epsc')

%% Comparison between output modes

figure('Name', 'Comparison of output modes', 'Position', [0 0 1700 700])
plot_ellipse(epsilon_oseen,theta_oseen,lcp_oseen,X(row_oseen(1),col_oseen),Y(row_oseen,col_oseen(1)),'blue','+',true)
%plot_ellipse(epsilon_topf,theta_topf,lcp_topf,X(row_topf(1),col_topf),Y(row_topf,col_topf(1)),'green','+',true)
plot_ellipse(epsilon,theta,lcp,X(row(1),col),Y(row,col(1)),'red','+',true)
title('Comparison of output modes')
xlabel('Detuning')
ylabel('Gain')

% Need some custom legend here
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'k-','LineWidth',2);
h(2) = plot(NaN,NaN,'k:','LineWidth',2);
h(3) = plot(NaN,NaN,'Color', 'blue','LineWidth',2);
%h(4) = plot(NaN,NaN,'Color', 'green','LineWidth',2);
h(4) = plot(NaN,NaN,'Color', 'red','LineWidth',2);
legend(h, 'RCP', 'LCP', 'Oseen', 'Dev. CWT')

daspect([1 1 1])
savefig(gcf, 'simple_cavity', 'Comparison modes', 'png')


