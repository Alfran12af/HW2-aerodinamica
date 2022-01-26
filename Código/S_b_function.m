function [S_b] = S_b_function(inc, aoa, r)
delta = inc+aoa;
x = sqrt(r^2/(r^2+(tan((pi/2)-delta))^2));
y = tan((pi/2)-delta)*x;

dist = sqrt(x^2+y^2);

S_b = pi*dist*r;
end