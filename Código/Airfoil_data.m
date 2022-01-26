clc
clear all
close all

load('airfoil_data_workspace.mat') % open imported airfoil (CH60) data from csv file


%% Extract csv data
A = table2array(mh60pol1);
Re = [100000; 200000; 400000; 600000; 1000000];
C_T = 0.28; % wing tip chord
C_R = 1.55; % wing root chord


%% Cl as a function of alpha
for i = 1:5 % compute every reynolds number
    cl_matrix(i,:) = A(i,2:8:end); % create cl matrix
    alpha_matrix(i,:) = A(i,1:8:end); %create alpha matrix
end


figure
hold on
title('$C_l$ as a function of $\alpha$', 'Interpreter','latex', FontSize=16)
plot(alpha_matrix(1,:), cl_matrix(1,:), '-.') % reynolds out of service range
plot(alpha_matrix(2,:), cl_matrix(2,:))
plot(alpha_matrix(3,:), cl_matrix(3,:))
plot(alpha_matrix(4,:), cl_matrix(4,:))
plot(alpha_matrix(5,:), cl_matrix(5,:))
grid on
grid minor
xlabel('$\alpha$ (deg.)', 'Interpreter','latex', FontSize=15)
ylabel('$C_l$', 'Interpreter','latex', FontSize=15)
legend('$Re=100000$', '$Re=200000$', '$Re=400000$', '$Re=600000$', '$Re=1000000$', 'location', 'southeast', 'Interpreter', 'latex', FontSize=12)
hold off


% linear regression for Re=600000
lin_reg_T = polyfit(A(4,1:8:176), A(4,2:8:176), 1); % wing tip regression line
yfit_T= polyval(lin_reg_T, A(4,1:8:176));
Cl_alpha_T = lin_reg_T(1,1);
alpha_l0_T = -lin_reg_T(1,2)/lin_reg_T(1,1);

% and Re=1000000
lin_reg_R = polyfit(A(5,1:8:176), A(5,2:8:176), 1); % wing root regression line
yfit_R= polyval(lin_reg_R, A(5,1:8:176));
Cl_alpha_R = lin_reg_R(1,1);
alpha_l0_R = -lin_reg_R(1,2)/lin_reg_R(1,1);

% for Re=600000

figure
hold on
title('$C_l$ as a function of $\alpha$ and regression line for Re=600000', 'Interpreter','latex', FontSize=16)
plot(alpha_matrix(4,:), cl_matrix(4,:), LineWidth=1.2)
plot(A(4,1:8:176), yfit_T, LineWidth=1.2)
plot(alpha_l0_T, 0, 'r*')
grid on
grid minor
xlabel('$\alpha$ (deg.)', 'Interpreter','latex', FontSize=15)
ylabel('$C_l$', 'Interpreter','latex', FontSize=15)
legend('$C_l$ curve', 'Regression line', 'Interpreter', 'latex', 'location', 'northwest', FontSize=12)
hold off

% for Re=1000000
figure
hold on
title('$C_l$ as a function of $\alpha$ and regression line for Re=1000000', 'Interpreter','latex', FontSize=16)
plot(alpha_matrix(5,:), cl_matrix(5,:), LineWidth=1.2)
plot(A(5,1:8:176), yfit_R, LineWidth=1.2)
plot(alpha_l0_R, 0, 'r*')
grid on
grid minor
xlabel('$\alpha$ (deg.)', 'Interpreter','latex', FontSize=15)
ylabel('$C_l$', 'Interpreter','latex', FontSize=15)
legend('$C_l$ curve', 'Regression line', 'Interpreter', 'latex', 'location', 'northwest', FontSize=12)
hold off


%% Cd as a function of Cl
for i = 1:5 % compute every reynolds number
    cd_matrix(i,:) = A(i,3:8:end); % create cd matrix
end


figure
hold on
title('$C_d$ as a function of $C_l$', 'Interpreter','latex', FontSize=16)
plot(cl_matrix(1,:), cd_matrix(1,:), '-.') % reynolds out of service range
plot(cl_matrix(2,:), cd_matrix(2,:))
plot(cl_matrix(3,:), cd_matrix(3,:))
plot(cl_matrix(4,:), cd_matrix(4,:))
plot(cl_matrix(5,:), cd_matrix(5,:))
grid on
grid minor
xlabel('$C_l$', 'Interpreter','latex', FontSize=15)
ylabel('$C_d$', 'Interpreter','latex', FontSize=15)
legend('$Re=100000$', '$Re=200000$', '$Re=400000$', '$Re=600000$', '$Re=1000000$','location', 'northwest', 'Interpreter', 'latex', FontSize=12)
hold off


% curva regresión para el cd
% Re=600000
polin_reg_T = polyfit(A(4,2:8:176), A(4,3:8:176), 2); % wing tip
pfit_T= polyval(polin_reg_T, A(4,2:8:176));
CD0_T = polin_reg_T(1,3);

figure
hold on
title('$C_d$ as a function of $C_l$', 'Interpreter','latex', FontSize=16)
plot(cl_matrix(4,:), cd_matrix(4,:), LineWidth=1.2)
plot(A(4,2:8:176), pfit_T, LineWidth=1.2)
grid on
grid minor
xlabel('$C_l$', 'Interpreter','latex', FontSize=15)
ylabel('$C_d$', 'Interpreter','latex', FontSize=15)
legend('$Re=600000$','location', 'northwest', 'Interpreter', 'latex', FontSize=12)
hold off

% Re=1000000
polin_reg_R = polyfit(A(5,2:8:176), A(5,3:8:176), 2); % wing root
pfit_R= polyval(polin_reg_R, A(5,2:8:176));
CD0_R = polin_reg_R(1,3);

figure
hold on
title('$C_d$ as a function of $C_l$', 'Interpreter','latex', FontSize=16)
plot(cl_matrix(5,:), cd_matrix(5,:), LineWidth=1.2)
plot(A(5,2:8:176), pfit_R, LineWidth=1.2)
grid on
grid minor
xlabel('$C_l$', 'Interpreter','latex', FontSize=15)
ylabel('$C_d$', 'Interpreter','latex', FontSize=15)
legend('$Re=1000000$','location', 'northwest', 'Interpreter', 'latex', FontSize=12)
hold off


%% Moment coefficient
% the given moment coefficient is refered to the aerodynamic center

for i = 1:5 % compute every reynolds number
    cm_ac_matrix(i,:) = A(i,4:8:end); % create cm matrix
end

% compute the moment coefficient from ac for every Reynolds number
cm_ac_matrix_T = cm_ac_matrix(4,:); % tip free moment coefficient
cm_ac_matrix_R = cm_ac_matrix(5,:); % root free moment coefficient

% compute the mean moment coefficient from ac for the wing tip
CM0_T = mean(cm_ac_matrix_T(1,1:22));

% compute the mean moment coefficient from ac for the wing root
CM0_R = mean(cm_ac_matrix_R(1,1:22));

% plot CM0 for the tip and root Re as a function of alpha
figure
hold on
title('$C_m$ as a function of $\alpha$ for wint root and tip', 'Interpreter','latex', FontSize=16)
plot(alpha_matrix(5,1:22), cm_ac_matrix_T(1,1:22))
plot(alpha_matrix(5,1:22), cm_ac_matrix_R(1,1:22))
grid on
grid minor
xlabel('$\alpha$', 'Interpreter','latex', FontSize=15)
ylabel('$C_{m_{ac}}$', 'Interpreter','latex', FontSize=15)
legend('$Re=600000$', '$Re=1000000$', 'location', 'northeast', 'Interpreter', 'latex', FontSize=12)
ylim([-0.8 0.8])
hold off
