

alfa = -10:0.1:10;
U_inf = 100/3.6;
rho_inf = 1.225;
l_F = 2;

fun = @(x) 0.6^2*(1-x.^2);

M_F = pi/4*k*rho_inf*U_inf*U_inf*alfa*integral(fun,0,l_F);


ld = 2/0.6;
Sb = 

C_D0 = Cf*(1+60/ld.^3+0.0025*ld)Ss/Sb;