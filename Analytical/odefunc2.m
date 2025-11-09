function dvar = odefunc2(~,var,miu)
L = 1;
r = 0.002;
m = 0.402;
J = m*L*L/12;
g = 9.81;

th = var(1);
dth = var(3);

A = J + (L*m*(L*cos(th) - 2*r*sin(th))*(cos(th) - miu*sin(th)))/4 - (m*r*(L*cos(th) - 2*r*sin(th))*(sin(th) + miu*cos(th)))/2;
B = (m*(L*sin(th) + 2*r*cos(th))*(2*r*sin(th) - L*cos(th) + L*miu*sin(th) + 2*miu*r*cos(th)))/4;
C = (g*m*(2*r*sin(th) - L*cos(th) + L*miu*sin(th) + 2*miu*r*cos(th)))/2;

d2th = (C-B*dth^2)/A;

dvar(1,1) = var(3);
dvar(2,1) = var(4);
dvar(3,1) = d2th;
dvar(4,1) = miu*(d2th*r*sin(th) - g + (L*dth^2*sin(th))/2 + dth^2*r*cos(th) - (L*d2th*cos(th))/2);