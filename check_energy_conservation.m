% Check energy conservation of the CWT.

close all;

%% Data for the medium
Lp = 300e-9; % m
epsilona = 3.2;
epsilonb = 2.9;
L = 20*Lp; % m
psi = 0%pi/4

n1 = sqrt((epsilona+epsilonb)/2);
n2 = n1;

subtitle = sprintf('$L_p=%.3g\\mu$m, $\\epsilon_a = %.3g+i%.3g$, $\\epsilon_b = %.3g+i%.3g$,\n $L=%.3g\\mu$m $n_1=%.3g$, $n_2=%.3g$', Lp*1e6, real(epsilona), imag(epsilona), real(epsilonb), imag(epsilonb), L*1e6, n1, n2);

save_figures = true;

%% Intermediary variables
p = 2*pi/Lp;
lambda_bragg = real(sqrt((epsilona+epsilonb)/2))*Lp;
delta_lambda = real(sqrt((epsilona+epsilonb)/2)-sqrt(epsilonb))*Lp;
lambda = linspace(lambda_bragg-3*delta_lambda,lambda_bragg+3*delta_lambda,1000);
k_0 = 2*pi./lambda;

%% Reflections and transmissions for CWT method

R_cwt = zeros(length(lambda), 4);
T_cwt = zeros(length(lambda), 4);

for j=1:length(lambda)
    M_interface_1 = interface_chiral_to_isotrope(p, epsilona, epsilonb, psi, 0, k_0(j), n1)^-1;
    M_medium = cwt(p, epsilona, epsilonb, psi, L, k_0(j));
    M_interface_2 = interface_chiral_to_isotrope(p, epsilona, epsilonb, psi, L, k_0(j), n2);
    M = M_interface_2*M_medium*M_interface_1;
    R = reflection_cwt(M);
    T = transmission_cwt(M);
    R_cwt(j,:) = reshape(R, [], 4);
    T_cwt(j,:) = reshape(T, [], 4);

    M_oseen = oseen(p, epsilona, epsilonb, psi, L, k_0(j));
    R = reflection_oseen(M_oseen, n1, n2);
    T = transmission_oseen(M_oseen, n1, n2);
    R_oseen(j,:) = reshape(R, [], 4);
    T_oseen(j,:) = reshape(T, [], 4);
end

r_LL = R_cwt(:,1);
r_LR = R_cwt(:,2);
r_RL = R_cwt(:,3);
r_RR = R_cwt(:,4);
t_LL = T_cwt(:,1);
t_LR = T_cwt(:,2);
t_RL = T_cwt(:,3);
t_RR = T_cwt(:,4);

conservation_L = abs(r_LL).^2 + abs(r_RL).^2 + (n2/n1)*(abs(t_LL).^2+abs(t_RL).^2) - 1;
conservation_R = abs(r_LR).^2 + abs(r_RR).^2 + (n2/n1)*(abs(t_LR).^2+abs(t_RR).^2) - 1;
redistribution = r_LL.*conj(r_LR) + r_RL.*conj(r_RR) + (n2/n1)*(t_LL.*conj(t_LR)+t_RL.*conj(t_RR));

r_LL_oseen = R_oseen(:,1);
r_LR_oseen = R_oseen(:,2);
r_RL_oseen = R_oseen(:,3);
r_RR_oseen = R_oseen(:,4);
t_LL_oseen = T_oseen(:,1);
t_LR_oseen = T_oseen(:,2);
t_RL_oseen = T_oseen(:,3);
t_RR_oseen = T_oseen(:,4);

conservation_L_oseen = abs(r_LL_oseen).^2 + abs(r_RL_oseen).^2 + (n2/n1)*(abs(t_LL_oseen).^2+abs(t_RL_oseen).^2) - 1;
conservation_R_oseen = abs(r_LR_oseen).^2 + abs(r_RR_oseen).^2 + (n2/n1)*(abs(t_LR_oseen).^2+abs(t_RR_oseen).^2) - 1;
redistribution_oseen = r_LL_oseen.*conj(r_LR_oseen) + r_RL_oseen.*conj(r_RR_oseen) + (n2/n1)*(t_LL_oseen.*conj(t_LR_oseen)+t_RL_oseen.*conj(t_RR_oseen));
%% Plot CWT method
f_cwt = figure('Name', 'CWT');
hold on
plot(lambda*1e9, conservation_L, '--', 'DisplayName','Conservation LH', 'LineWidth', 2)
plot(lambda*1e9, conservation_R, '-', 'DisplayName','Conservation RH', 'LineWidth', 2)
plot(lambda*1e9, abs(redistribution), '-.', 'DisplayName','$|\mathrm{Redistribution}|$', 'LineWidth', 2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (m)', 'interpreter', 'latex')
ylabel('Conservation', 'interpreter', 'latex')
hold off

%% Plot Oseen method
f_oseen = figure('Name', 'Oseen');
hold on
plot(lambda*1e9, conservation_L_oseen, '--', 'DisplayName','Conservation LH', 'LineWidth', 2)
plot(lambda*1e9, conservation_R_oseen, '-', 'DisplayName','Conservation RH', 'LineWidth', 2)
plot(lambda*1e9, abs(redistribution_oseen), '-.', 'DisplayName','$|\mathrm{Redistribution}|$', 'LineWidth', 2)
legend('interpreter', 'latex')
title(subtitle, 'interpreter', 'latex')
xlabel('$\lambda$ (m)', 'interpreter', 'latex')
ylabel('Conservation', 'interpreter', 'latex')
hold off

%% Save figures

if save_figures
savefig(f_oseen, 'energy_conservation', 'oseen', 'eps', 'epsc')
savefig(f_cwt, 'energy_conservation', 'cwt', 'eps', 'epsc')
end
