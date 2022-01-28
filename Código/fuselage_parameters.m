function [C_MF_adim, C_MF_2_adim, C_D0_adim] = fuselage_parameters(rho_inf, U_inf, k, l_f, r, S, mac, span, ALPHA, inc)

q_inf = 1/2*rho_inf*U_inf^2;% dynamic pressure

% elipsoid area
p = @(x) (2*pi*sqrt(0.3^2*(1-x.^2)));
S_f = integral(p, -l_f, l_f);
S_s = S_f;

syms alpha_inc % declare alpha as a variable

%% Pitching moment
b_F = @(x) (r^2*(1-x.^2));
I1 = integral(b_F, -l_f, l_f);

Sol_M_F = pi/2*k*q_inf*alpha_inc*4*I1;

M_F = pi/2*k*q_inf*alpha_inc*4*I1;

alpha_vec  = [-3*pi/180:0.5*pi/180:10*pi/180];
M_F_vec = pi/2*k*q_inf.*alpha_vec(1,:)*4*I1;
alpha_vec_deg = [-3:0.5:10];
C_MF = M_F_vec/(q_inf*S_f*2*l_f);

% calculamos el CM_F adimensional respecto el ala, con
span = 20;
S_alar = S*span^2;
mac_alar = mac*span;
C_MF_adim = M_F_vec/(q_inf*S_alar*mac_alar);

% calculamos el coeficiente de momento libre para el fuselaje. Usamos los vectroes alpha para el apartado 2 
M_F_vec_2 = pi/2*k*q_inf.*(ALPHA(1,:)*pi/180)*4*I1;
C_MF_2_adim = M_F_vec_2/(q_inf*S_alar*mac_alar); % subindex 2 indicates that is for computing coeff parameters

% set color for plots
str1 = '#77AC30';
green = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#EDB120';
yellow = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
str3 = '#7E2F8E';
lilac = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;
str4 = '#D95319';
orange = sscanf(str4(2:end),'%2x%2x%2x',[1 3])/255;


figure
hold on
title('Fuselage pitching moment coefficient as a function of $\alpha$', 'interpreter', 'latex', 'FontSize',16)
plot(alpha_vec_deg, C_MF_adim, 'Color', green, 'LineWidth',1.2)
xlabel('alpha ', 'Interpreter','latex', 'FontSize',15)
ylabel('Fuselage pitching moment coefficient ($C_{MF}$)', 'Interpreter','latex', 'FontSize',15)
grid on
grid minor
axis padded
hold off

%% Parasitic drag
aoa=0*pi/180;

[S_b] = S_b_function(inc, aoa, r);

c_f = 0.00361985; % friction coefficient
C_Db = 0; % cause we have a closed fuselage
C_D0 = c_f*(1+60/(2*l_f/(2*r))^3+0.0025*(2*l_f/(2*r)))*S_s/S_b + C_Db; % parasitic drag coefficient
D = q_inf*C_D0*S_b; % total parasitic drag [N]

% adimensionalizamos el coeficiente de drag respecto del ala
C_D0_adim = (c_f*(1+60/((2*l_f/(2*r))^3)+0.0025*(2*l_f/(2*r)))*S_s/S_b + C_Db)*(S_b/(S_alar));

% compute the CD0 as a function of alpha, we define it as C_Dp
for i=1:length(alpha_vec)
    [S_b_alpha] = S_b_function(inc, alpha_vec(1,i), r); % S_b as a function of alpha
    
    C_Dp(i) = (c_f*(1+60/(2*l_f/(2*r))^3+0.0025*(2*l_f/(2*r))).*S_s/S_b_alpha + C_Db)*(S_b/(S_alar)); % parasitic drag vector adimensional with the wing
end

figure
hold on
title('Fuselage parasitic drag coefficient as a function of $\alpha$', 'interpreter', 'latex', 'FontSize',16)
plot(alpha_vec_deg, C_Dp(1,:), 'Color', orange, 'LineWidth',1.2)
xlabel('alpha [deg.]', 'Interpreter','latex', 'FontSize',15)
ylabel('Fuselage parasitic drag coefficient ($C_{D0}$)', 'Interpreter','latex', 'FontSize',15)
grid on
grid minor
axis padded
hold off

end