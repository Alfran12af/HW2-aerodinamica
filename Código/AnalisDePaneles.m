% ///////////////////////////////////////////////////////////
% / Sript used to calculate to make the analysis of Panels  /        /
% /                                                         /
% / GRETA                                                   /
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

Paneles = 150;
RangPaneles = (1:1:Paneles);
ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ]
alpha_l0_Mat = zeros(1,Paneles);

for i = RangPaneles
    [CM_cg_matrix, force_coeff] = LLWing_function(0,i,1);
    Cl_matrix = force_coeff(7,:);
    lin_reg = polyfit(ALPHA(1,:), Cl_matrix(1,:), 1);
    alpha_l0 = -lin_reg(1,2)/lin_reg(1,1);
    alpha_l0_Mat(1,i) = alpha_l0;
end

figure
hold on
title('Análisis del angulo de sustentación nula en función de los paneles', 'Interpreter', 'latex', 'FontSize',16)
plot(RangPaneles, alpha_l0_Mat, 'LineWidth',1.2)
xlabel ('Numero de Paneles')
ylabel ('Angulo de sustentacion nula')
ylim([0.8,1.5])
grid on
grid minor    
    
    
    
