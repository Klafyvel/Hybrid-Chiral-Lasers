% Simulates a medium composed of two right handed slab, but one is twisted. Use Oseen and CWT method and compare them.
close all;
%% Data for the medium
Lp = 300e-9; % m
epsilona = 3.2%+0.02i;
epsilonb = 2.9%+0.02i;
L = 10*Lp; % m
twist = 0;%pi/3;

n1 = 1;%sqrt((epsilona+epsilonb)/2);
n2 = 1;%n1;

subtitle = sprintf('$L_p=%.3g\\mu$m, $\\epsilon_a = %.3g+i%.3g$, $\\epsilon_b = %.3g+i%.3g$,\n $L=%.3g\\mu$m $n_1=%.3g$, $n_2=%.3g$, $\\delta\\psi=%.3g^{\\circ}$', Lp*1e6, real(epsilona), imag(epsilona), real(epsilonb), imag(epsilonb), L*1e6, n1, n2, twist*180/pi);

lhm = false;

%% Useful intermediary variables
margin=0.01;
lambda_bragg = real(sqrt((epsilona+epsilonb)/2))*Lp;
delta_lambda = 2*real(sqrt((epsilona+epsilonb)/2)-sqrt(epsilonb))*Lp;
lambda = [linspace(lambda_bragg-3*delta_lambda,lambda_bragg-margin*delta_lambda,2000) linspace(lambda_bragg-margin*delta_lambda,lambda_bragg+margin*delta_lambda,1000) linspace(lambda_bragg+margin*delta_lambda,lambda_bragg+3*delta_lambda,2000)];
if ~lhm
    p = 2*pi/Lp;
else
    p = -2*pi/Lp;
end
k0 = 2.*pi./lambda;
comp = zeros(1, length(lambda));

%% Reflectivitys and transmissions for oseen method

R_oseen = zeros(length(lambda), 4);
T_oseen = zeros(length(lambda), 4);
R_cwt = zeros(length(lambda), 4);
T_cwt = zeros(length(lambda), 4);

for j=1:length(lambda)
    M_oseen1 = oseen(p, epsilona, epsilonb, 0, L, k0(j));
    M_oseen2 = oseen(p, epsilona, epsilonb, twist, L, k0(j));
    R = reflection_oseen(M_oseen2*M_oseen1, n1, n2);
    T = transmission_oseen(M_oseen2*M_oseen1, n1, n2);
    R_oseen(j,:) = reshape(R, [], 4);
    T_oseen(j,:) = reshape(T, [], 4);

    M_interface_1 = interface_chiral_to_isotrope(p, epsilona, epsilonb, 0, 0, k0(j), n1, lhm)^-1;
    M_medium_1 = cwt(p, epsilona, epsilonb, 0, L, k0(j), lhm);
    M_interface_2 = interface_chiral_same_handedness(p, epsilona, epsilonb, 0, p, epsilona, epsilonb, twist, L, k0(j), lhm);
    M_medium_2 = cwt(p, epsilona, epsilonb, twist, L, k0(j), lhm);
    M_interface_3 = interface_chiral_to_isotrope(p, epsilona, epsilonb, twist, L, k0(j), n2, lhm);
    M = M_interface_3*M_medium_2*M_interface_2*M_medium_1*M_interface_1;
    R = reflection_cwt(M);
    T = transmission_cwt(M);
    R_cwt(j,:) = reshape(R, [], 4);
    T_cwt(j,:) = reshape(T, [], 4);
end
% Calculate intensity reflections
R_oseen = abs(R_oseen) .^ 2;
T_oseen = n2/n1*abs(T_oseen) .^ 2;

% Calculate intensity reflections
R_cwt = abs(R_cwt) .^ 2;
T_cwt = n2/n1*abs(T_cwt) .^ 2;

% Calculate comparison
R_comp = R_oseen - R_cwt;
T_comp = T_oseen - T_cwt;

%% BEYOND THIS POINT STANDS THE UNHOLY LAND OF PLOTTING. ADVANCE AT YOUR OWN RISKS.


%% Plot Oseen method
figure('Name','Reflectivity, Exact theory')
hold on
plot(lambda.*1e9, R_oseen(:,1), '--', 'DisplayName','$R_{LL}$','LineWidth',2)
plot(lambda.*1e9, R_oseen(:,2), ':', 'DisplayName','$R_{LR}$','LineWidth',2)
plot(lambda.*1e9, R_oseen(:,3), '-.', 'DisplayName','$R_{RL}$','LineWidth',2)
plot(lambda.*1e9, R_oseen(:,4), 'DisplayName','$R_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'oseen_reflection', 'png')
savefig(gcf, 'narrow_band', 'oseen_reflection', 'eps', 'epsc')
hold off
figure('Name', 'Transmissivity, Exact theory')
hold on
plot(lambda.*1e9, T_oseen(:,1), '--', 'DisplayName','$T_{LL}$','LineWidth',2)
plot(lambda.*1e9, T_oseen(:,2), ':', 'DisplayName','$T_{LR}$','LineWidth',2)
plot(lambda.*1e9, T_oseen(:,3), '-.', 'DisplayName','$T_{RL}$','LineWidth',2)
plot(lambda.*1e9, T_oseen(:,4), 'DisplayName','$T_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmissivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'oseen_transmission', 'png')
savefig(gcf, 'narrow_band', 'oseen_transmission', 'eps', 'epsc')
hold off



%% Plot CWT method
figure('Name', 'Reflectivity, CWT method')
hold on
plot(lambda.*1e9, R_cwt(:,1), '--', 'DisplayName','$R_{LL}$','LineWidth',2)
plot(lambda.*1e9, R_cwt(:,2), ':', 'DisplayName','$R_{LR}$','LineWidth',2)
plot(lambda.*1e9, R_cwt(:,3), '-.', 'DisplayName','$R_{RL}$','LineWidth',2)
plot(lambda.*1e9, R_cwt(:,4), 'DisplayName','$R_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'cwt_reflection', 'png')
savefig(gcf, 'narrow_band', 'cwt_reflection', 'eps', 'epsc')
hold off
figure('Name', 'Transmissivity, CWT method')
hold on
plot(lambda.*1e9, T_cwt(:,1), '--', 'DisplayName','$T_{LL}$','LineWidth',2)
plot(lambda.*1e9, T_cwt(:,2), ':', 'DisplayName','$T_{LR}$','LineWidth',2)
plot(lambda.*1e9, T_cwt(:,3), '-.', 'DisplayName','$T_{RL}$','LineWidth',2)
plot(lambda.*1e9, T_cwt(:,4), 'DisplayName','$T_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmissivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'cwt_transmission', 'png')
savefig(gcf, 'narrow_band', 'cwt_transmission', 'eps', 'epsc')
hold off


%% Plot comparison
figure('Name', 'Reflectivity, difference between CWT and Oseen')
hold on
plot(lambda.*1e9, R_comp(:,1), '--', 'DisplayName','$R_{LL}$','LineWidth',2)
plot(lambda.*1e9, R_comp(:,2), ':', 'DisplayName','$R_{LR}$','LineWidth',2)
plot(lambda.*1e9, R_comp(:,3), '-.', 'DisplayName','$R_{RL}$','LineWidth',2)
plot(lambda.*1e9, R_comp(:,4), 'DisplayName','$R_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'comparison_reflection', 'png')
savefig(gcf, 'narrow_band', 'comparison_reflection', 'eps', 'epsc')
hold off

figure('Name', 'Transmissivity, difference between CWT and Oseen')
hold on
plot(lambda.*1e9, T_comp(:,1), '--', 'DisplayName','$T_{LL}$','LineWidth',2)
plot(lambda.*1e9, T_comp(:,2), ':', 'DisplayName','$T_{LR}$','LineWidth',2)
plot(lambda.*1e9, T_comp(:,3), '-.', 'DisplayName','$T_{RL}$','LineWidth',2)
plot(lambda.*1e9, T_comp(:,4), 'DisplayName','$T_{RR}$','LineWidth',2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmissivity', 'interpreter', 'latex')
savefig(gcf, 'narrow_band', 'comparison_transmission', 'png')
savefig(gcf, 'narrow_band', 'comparison_transmission', 'eps', 'epsc')
hold off

