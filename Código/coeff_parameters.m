% ///////////////////////////////////////////////////////////
% / Using cl_local, function that calculates the            /
% / cl_alpha and alpha_l0 and plots the regression line     /
% / in both cases                                           /
% ///////////////////////////////////////////////////////////

function [alpha_l0, cl_alpha, mom_slope, x_ac, CM0, CLa, CD0, k1, mean_CM_cg, CM_cg_matrix, CL_trim, CL_stall, v_stall, margin_stab] = coeff_parameters(ALPHA, force_coeff, cl_local, S, span, mac, C_MF_2_adim, C_D0_adim, rho_inf, BreakCl)

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

% take the lift coefficient of the wing for each aoa
Cl_matrix = force_coeff(7,:);

% using regression line, compute cl_alpha and alpha_l0
lin_reg = polyfit(ALPHA(1,:), Cl_matrix(1,:), 1); % wing regression line
%yfit = polyval(lin_reg, ALPHA(1,:));
cl_alpha = lin_reg(1,1);
alpha_l0 = -lin_reg(1,2)/lin_reg(1,1);

if BreakCl==1
    cl_alpha = 0; 
    mom_slope = 0; 
    x_ac = 0; 
    CM0 = 0; 
    CLa = 0; 
    CD0 = 0; 
    k1 = 0; 
    mean_CM_cg = 0; 
    CM_cg_matrix = 0; 
    CL_trim = 0; 
    CL_stall = 0; 
    v_stall = 0; 
    margin_stab = 0;
    return
end 
figure
hold on
title('Regression line for $C_l$ as a function of $\alpha$', 'Interpreter','latex', 'FontSize',16)
plot(ALPHA(1,:), Cl_matrix(1,:), 'Color', yellow, 'LineWidth',1.2)
plot(ALPHA(1,:), Cl_matrix(1,:), 'bo', 'MarkerSize',8)
% plot(ALPHA, yfit, LineWidth=1.2)
plot(alpha_l0, 0, 'r*')
grid on
grid minor
xlabel('$\alpha$ (deg.)', 'Interpreter','latex', 'FontSize',15)
ylabel('$C_l$', 'Interpreter','latex', 'FontSize',15)
hold off


%% compute the moment coefficient of the wing
CMF = C_MF_2_adim;
CMw_LE = force_coeff(5,:);
CM_tot = CMw_LE;

% using regression line, compute the slope of moment as a function of alpha
mom_lin_reg = polyfit(force_coeff(7,:), force_coeff(5,:), 1); % wing regression line
%mom_yfit = polyval(mom_lin_reg, Cl_matrix);
mom_slope = mom_lin_reg(1,1);
% as we compute x_ac from the lift and moment coefficient, it is adimensionalized with the mac
x_ac  = -mom_lin_reg(1,1)*mac*span; % compute the real x_ac from the leadin edge
CM0 = force_coeff(5,:)+force_coeff(7,:)*x_ac;


figure
hold on
title('Regression line for $C_{M_{LE}}$ as a function of $C_L$', 'Interpreter', 'latex', 'FontSize',16)
plot(Cl_matrix(1,:), CM_tot(:), 'Color', red, 'LineWidth',1.2)
plot(Cl_matrix(1,:), CM_tot(:), 'bo', 'MarkerSize',8)
%plot(0, CM0, 'r*')
grid on
grid minor
xlabel('$C_L$', 'Interpreter','latex', 'FontSize',15)
ylabel('$C_{M_{LE}}$', 'Interpreter','latex', 'FontSize',15)
hold off


%% Basic and aditional lift coefficient distributions
% lo haremos mediante la resolución de dos ecuaciones, dados dos angulos de ataque y dos cl
% los angulos de ataque seleccionados son los correspondientes a -4º y 4º

CL1 = force_coeff(7,3); % wing cl corresponding to -4º
CL2 = force_coeff(7,5); % wing cl corresponding to 4º

%Sacamos el cl local para -4º y 4º
Cl1 = cl_local(:,3); % local cl for -4º
Cl2 = cl_local(:,5); % local cl for 4º

for i = 1:length(cl_local)
    CLa(i,1) = (Cl1(i,1)-Cl2(i,1))/(CL1-CL2);
end
for i = 1:length(cl_local)
    CLb(i,1) = Cl1(i,1)-CLa(i,1)*CL1;
end

% create y matrix for plot
y = linspace(-10, 10, length(cl_local));

% cl distribution
found = 0;
cl_max_T = 1.228;
for i = 1:length(CLa)
    cl_y(i,1) = (cl_max_T-CLb(i,1))/CLa(i,1); % cl max of the wing tip 
end
cl_max = min(cl_y);

CL_stall = -1;
while found ~= 1
    CL_stall = CL_stall + 0.01;

    cl_dist = CLa*CL_stall+CLb;
    if (max(cl_dist)>=cl_max)
        found = 1;
    end  
end

cl_max_vec = ones(1,length(cl_local)).*cl_max;

% compute the stall speed
W_S = 20; % weight surface ratio [kg/m^2]
g = 9.81; % gravity acceleration [m/s^2]
v_stall = sqrt(2*W_S*g/(rho_inf*cl_max));


figure
hold on
title('Cl básico y Cl adicional a lo largo de la envergadura', 'Interpreter','latex', 'FontSize',16)
plot(y(1,:), CLa(:,1), 'Color', lilac, 'LineWidth',1.2)
plot(y(1,:), CLb(:,1), 'b', 'LineWidth',1.2)
plot(y(1,:), cl_max_vec, 'Color', orange, 'LineWidth',1.2)
plot(y(1,:), cl_dist, 'Color', green, 'LineWidth',1.2)
grid on
grid minor
xlabel('Distancia y desde la raiz', 'Interpreter', 'latex', 'FontSize',15)
ylabel('$C_{l}$', 'Interpreter','latex', 'FontSize',15)
legend('$C_{la}$', '$C_{lb}$', '$C_{lmax}$', '$ distribucion C_{l} $', 'Interpreter','latex', 'FontSize',12')
hold off

%% Calculate CD0 of all wing+fuselage
CD0_w = force_coeff(11,:); % drag coefficient of the wing (profile+induced)
CD0_f = C_D0_adim; % drag coefficient of the fuselage
CD0_matrix = CD0_w+CD0_f;

polin_reg = polyfit(Cl_matrix.^2, CD0_matrix, 1); % polar regresion for CD0 of the wing and fuselage
pfit_CD0 = polyval(polin_reg, linspace(-1, 1, 80));

CD0 = polin_reg(1,2);
k1 = polin_reg(1,1);
polin_reg2 = polyfit(Cl_matrix,CD0_matrix,2);
polin_reg2(2) = 0;
pfit_CD02 = polyval(polin_reg2, linspace(-1, 1, 80));

polin_reg3 = polyfit(Cl_matrix,CD0_w,2);
polin_reg3(2) = 0;
pfit_CD03 = polyval(polin_reg3, linspace(-1, 1, 80));

%k2 = polin_reg(1,2);

%Drag parabolic curve
  
disp(polin_reg3)
disp(polin_reg2)


figure
hold on
title('Distribución del Drag coeficient en funcion del Cl', 'Interpreter','latex', 'FontSize',16)
plot(linspace(-1, 1, 80), pfit_CD02, 'Color', orange, 'LineWidth',1.2)
plot(linspace(-1, 1, 80), pfit_CD03, 'Color', blue, 'LineWidth',1.2)
%plot(Cl_matrix(1,:), CD0_matrix(1,:), 'ko', 'MarkerSize',8)
%plot(Cl_matrix(1,:), CD0_w(1,:), 'ko', 'MarkerSize',8)
grid on
grid minor
legend('Fuselage Drag + Wing Drag','Wing Drag')
xlabel('$C_L$', 'Interpreter','latex', 'FontSize',15)
ylabel('$C_{D0}$', 'Interpreter','latex', 'FontSize',15)
hold off
xlim([-1 1])

%% CM about the center of mass as a function of lift coefficient
x_cg = 1.3/(span*mac);
% CM0 is the wing pitching moment, and CMF is the fuselage pitching moment, so to compute the CM_cg_matrix, we have to sum both moments to find the global pitching moment
CM_cg_matrix = CM0-Cl_matrix(1,:)*(x_ac-x_cg)+CMF;
mean_CM_cg = mean(CM_cg_matrix);

% stability margin
margin_stab = (x_ac-1.3)/1.3;

% compute the lift coefficient which has a 0 moment from the cg
mom_cg_lin_reg = polyfit(Cl_matrix(1,:), CM_cg_matrix(1,:), 1); % CM from xcg regression line
yfit= polyval( mom_cg_lin_reg, Cl_matrix(1,:));
CL_trim = -mom_cg_lin_reg(1,2)/mom_cg_lin_reg(1,1);

figure
hold on
title('Coeficiente de momento respecto el centro de gravedad en función del $C_L$', 'Interpreter','latex', 'FontSize',15)
%plot(Cl_matrix(1,:), CM_cg_matrix(1,:), 'Color', green, LineWidth=1.2)
%plot(Cl_matrix(1,:), CM_cg_matrix(1,:), 'o', MarkerSize=8)
plot(Cl_matrix(1,:), yfit, 'b', 'LineWidth',1.2)
plot(CL_trim, 0, 'r*')
grid on
grid minor
xlabel('$C_L$', 'Interpreter','latex', 'FontSize',15)
ylabel('$C_{M_{cg}}$', 'Interpreter','latex', 'FontSize',15)
hold off


%% L/D as a function of flight speed
v_matrix = 10:1:50; % speed matrix [m/s]

for i = 1:length(v_matrix)
    CL_mat(1,i) = 2*W_S*g/(rho_inf*v_matrix(1,i)^2);
    CD_mat(1,i) = k1*CL_mat(1,i)^2+CD0;
end
L_D = CL_mat./CD_mat;

figure
hold on
title('Lift-Frag ratio ($L/D$) as a function of flight velocity', 'Interpreter','latex', 'FontSize',16)
plot(v_matrix(1,:), L_D(1,:), 'b', 'LineWidth',1.2)
grid on
grid minor
xlabel('$v [m/s]$', 'Interpreter','latex', 'FontSize',15)
ylabel('$\frac{L}{D}$', 'Interpreter','latex', 'FontSize',15)
hold off

end