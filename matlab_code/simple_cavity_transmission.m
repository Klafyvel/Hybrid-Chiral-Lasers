%close all;
%clear all;
%% Data for the medium
Lp = 300e-9; % m
epsilona = 3.2+0.02i;
epsilonb = 2.9+0.02i;
epsilonc = 2.8+0.02i;
L = 20*Lp;
psi = 0; 
chi = 30*pi/180;

n1 = 1;%sqrt((epsilona+epsilonb)/2);
n2 = 2;%n1;

lhm = false;

save_figures = true;

%% Intermediary variables

epsilona = epsilona*epsilonc / (epsilonc*cos(chi)^2+epsilona*sin(chi)^2);

n_bar = (sqrt(epsilona)+sqrt(epsilonb))/2;%sqrt((epsilona+epsilonb)/2);
delta_n = sqrt(epsilona)-sqrt(epsilonb);

lambda_0_bragg = real(n_bar)*Lp;
delta_lambda_0_bragg = real(delta_n)*Lp;

lambda_0 = linspace(lambda_0_bragg-3*delta_lambda_0_bragg, lambda_0_bragg+3*delta_lambda_0_bragg, 1000);

if ~lhm
    p = 2*pi/Lp;
else
    p = -2*pi/Lp;
end

k_0 = (2*pi) ./ lambda_0;

%% Reflections and transmissions

R_oseen = zeros(length(lambda_0), 4);
T_oseen = zeros(length(lambda_0), 4);
R_cwt = zeros(length(lambda_0), 4);
T_cwt = zeros(length(lambda_0), 4);
comp = zeros(1, length(lambda_0));

for j=1:length(lambda_0)
    M_oseen = oseen(p, epsilona, epsilonb, psi, L, k_0(j));
    R = reflection_oseen(M_oseen, n1, n2);
    T = transmission_oseen(M_oseen, n1, n2);
    R_oseen(j,:) = reshape(R, [], 4);
    T_oseen(j,:) = reshape(T, [], 4);

    M_interface_1 = interface_chiral_to_isotrope(p, epsilona, epsilonb, psi, 0, k_0(j), n1, lhm)^-1;
    M_medium = cwt(p, epsilona, epsilonb, psi, L, k_0(j), lhm);
    M_interface_2 = interface_chiral_to_isotrope(p, epsilona, epsilonb, psi, L, k_0(j), n2, lhm);
    M = M_interface_2*M_medium*M_interface_1;
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

R_comp = R_cwt - R_oseen;
T_comp = T_cwt - T_oseen;
%% Plot
f1 = figure('Name', 'Reflection coefficients, Oseen method');
plot(lambda_0.*1e9, R_oseen(:,1), '--', lambda_0.*1e9, R_oseen(:,2), ':', lambda_0.*1e9, R_oseen(:,3), '-.', lambda_0.*1e9, R_oseen(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),0,1.1*max(R_oseen, [], 'all')])
legend('$R_{RR}$','$R_{LR}$','$R_{RL}$','$R_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')

f2 = figure('Name', 'Reflection coefficients, CWT method');
plot(lambda_0.*1e9, R_cwt(:,1), '--', lambda_0.*1e9, R_cwt(:,2), ':', lambda_0.*1e9, R_cwt(:,3), '-.', lambda_0.*1e9, R_cwt(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),0,1.1*max(R_cwt, [], 'all')])
legend('$R_{RR}$','$R_{LR}$','$R_{RL}$','$R_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')

f3 = figure('Name', 'Transmission coefficients, Oseen method');
plot(lambda_0.*1e9, T_oseen(:,1), '--', lambda_0.*1e9, T_oseen(:,2), ':', lambda_0.*1e9, T_oseen(:,3), '-.', lambda_0.*1e9, T_oseen(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),0,1.1*max(T_oseen, [], 'all')])
legend('$T_{RR}$','$T_{LR}$','$T_{RL}$','$T_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmittivity', 'interpreter', 'latex')

f4 = figure('Name', 'Transmission coefficients, CWT method');
plot(lambda_0.*1e9, T_cwt(:,1), '--', lambda_0.*1e9, T_cwt(:,2), ':', lambda_0.*1e9, T_cwt(:,3), '-.', lambda_0.*1e9, T_cwt(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),0,1.1*max(T_cwt, [], 'all')])
legend('$T_{RR}$','$T_{LR}$','$T_{RL}$','$T_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmittivity', 'interpreter', 'latex')

% Plot comparison
f5 = figure('Name','Reflection coefficients, difference between CWT and Oseen');
plot(lambda_0.*1e9, R_comp(:,1), '--', lambda_0.*1e9, R_comp(:,2), ':', lambda_0.*1e9, R_comp(:,3), '-.', lambda_0.*1e9, R_comp(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),1.1*min(R_comp, [], 'all'),1.1*max(R_comp, [], 'all')])
legend('$R_{RR}$','$R_{LR}$','$R_{RL}$','$R_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Reflectivity', 'interpreter', 'latex')

f6 = figure('Name','Transmission coefficients, difference between CWT and Oseen');
plot(lambda_0.*1e9, T_comp(:,1), '--', lambda_0.*1e9, T_comp(:,2), ':', lambda_0.*1e9, T_comp(:,3), '-.', lambda_0.*1e9, T_comp(:,4), 'LineWidth',2)
axis([min(lambda_0*1e9),max(lambda_0*1e9),1.1*min(T_comp, [], 'all'),1.1*max(T_comp, [], 'all')])
legend('$T_{RR}$','$T_{LR}$','$T_{RL}$','$T_{RR}$', 'interpreter', 'latex')
xlabel('$\lambda$ (nm)', 'interpreter', 'latex')
ylabel('Transmittivity', 'interpreter', 'latex')
%% Save figures
if save_figures
savefig(f1, 'simple_cavity', 'reflection_oseen', 'eps', 'epsc')
savefig(f2, 'simple_cavity', 'reflection_cwt', 'eps', 'epsc')
savefig(f3, 'simple_cavity', 'transmission_oseen', 'eps', 'epsc')
savefig(f4, 'simple_cavity', 'transmission_cwt', 'eps', 'epsc')
savefig(f5, 'simple_cavity', 'reflection_comp', 'eps', 'epsc')
savefig(f6, 'simple_cavity', 'transmission_comp', 'eps', 'epsc')
end
