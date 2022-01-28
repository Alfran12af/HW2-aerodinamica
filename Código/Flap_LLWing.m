% ///////////////////////////////////////////////////////////
% / Sript used to calculate the cmcg as a function          /
% / of the flap deflection                                  /
% / GRETA - Marc Lahoz                                      /
% ///////////////////////////////////////////////////////////

clear all
clc
close all


% set color for plots
str1 = '#77AC30';
green = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#EDB120';
yellow = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
str3 = '#7E2F8E';
lilac = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;
str4 = '#D95319';
orange = sscanf(str4(2:end),'%2x%2x%2x',[1 3])/255;
str5 = '#0072BD';
blue = sscanf(str5(2:end),'%2x%2x%2x',[1 3])/255;
str6 = '#A2142F';
red = sscanf(str6(2:end),'%2x%2x%2x',[1 3])/255;


DE_flap_matrix = (-20:5:20); % flap deflection (deg, positive:down)
N = 100;
BreakCl = 0;

for i = 1:length(DE_flap_matrix)
    [CM_cg_matrix, force_coeff] = LLWing_function(DE_flap_matrix(1,i),N, BreakCl);
    CM_cg(i,:) = CM_cg_matrix(1,:);
    CL_mat(i,:) = force_coeff(7,:);
end

k1 = 0.0323987330292590;
CD0 = 0.00744322797954547;

CL_max_range = sqrt(CD0/k1);
CL_min_speed = sqrt(3*CD0/k1);



figure
hold on
title('Center of mass moment coefficient ($C_{M_{cg}}$) as a function of flap deflection ($\eta$)', 'Interpreter','latex', 'FontSize',16)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(1,:),1),linspace(-2,2,30)), '--','LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(2,:),1),linspace(-2,2,30)), '--','LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(3,:),1),linspace(-2,2,30)), '--','LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(4,:),1),linspace(-2,2,30)), '--','LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(5,:),1),linspace(-2,2,30)), 'k', 'LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(6,:),1),linspace(-2,2,30)), '-.', 'LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(7,:),1),linspace(-2,2,30)), '-.', 'LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(8,:),1),linspace(-2,2,30)), '-.', 'LineWidth',1.2)
% plot(linspace(-2,2,30), polyval(polyfit(CL_mat(1,:),CM_cg(9,:),1),linspace(-2,2,30)), '-.', 'LineWidth',1.2)
plot(CL_mat(1,:), CM_cg(1,:), '--', 'LineWidth',1.2)
plot(CL_mat(2,:), CM_cg(2,:), '--', 'LineWidth',1.2)
plot(CL_mat(3,:), CM_cg(3,:), '--', 'LineWidth',1.2)
plot(CL_mat(4,:), CM_cg(4,:), '--', 'LineWidth',1.2)
plot(CL_mat(5,:), CM_cg(5,:), 'k', 'LineWidth',1.2)
plot(CL_mat(6,:), CM_cg(6,:), '-.', 'LineWidth',1.2)
plot(CL_mat(7,:), CM_cg(7,:), '-.', 'LineWidth',1.2)
plot(CL_mat(8,:), CM_cg(8,:), '-.', 'LineWidth',1.2)
plot(CL_mat(9,:), CM_cg(9,:), '-.', 'LineWidth',1.2)
grid on
grid minor
ylim([-0.7 0.7])
xlabel('$C_{L}$', 'Interpreter','latex', 'FontSize',15)
ylabel('$C_{M_{cg}}$', 'Interpreter','latex', 'FontSize',15)

% plot clmax lines
plot([CL_max_range CL_max_range], [-0.7 0.7], '--k', 'LineWidth',1.2);
plot([CL_min_speed CL_min_speed], [-0.7 0.7], ':k', 'LineWidth',1.2);
h1 = text(CL_max_range-0.05, 0.18, "Max. range", 'FontSize',12);
set(h1, 'Rotation', 90);
h2 = text(CL_min_speed-0.05, 0.18, "Min. speed", 'FontSize',12);
set(h2, 'Rotation', 90);
legend('$\eta = -20$', '$\eta = -15$', '$\eta = -10$', '$\eta = -5$', '$\eta = 0$', '$\eta = 5$', '$\eta = 10$', '$\eta = 15$', '$\eta = 20$', 'Max. range', 'Min. speed', 'Location', 'northeast', 'Interpreter','latex', 'FontSize',12)
hold off


