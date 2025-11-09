function dvar = odefunc1(~,var,miu)
L = 1;
r = 0.002;
m = 0.402;
J = m*L*L/12;
g = 9.81;
th = var(1);
dth = var(2);

A = (- m*L^2*miu^2 + m*L^2 - 4*m*cos(2*th)*L*miu*r - 4*m*sin(2*th)*L*r + 4*m*miu^2*r^2 + 4*J*miu^2 + 4*m*r^2 + 4*J)/(4*(miu^2 + 1));
B = (L*m*(L*miu - 2*r*cos(2*th) + 2*miu*r*sin(2*th)))/(2*(miu^2 + 1));
C = (g*m*(2*r*sin(th) - L*cos(th) + 2*L*miu*sin(th) + L*miu^2*cos(th) + 2*miu^2*r*sin(th)))/(2*(miu^2 + 1));

dvar(1,1) = var(2);
dvar(2,1) = (C-B*dth^2)/A;