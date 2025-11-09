% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

function dvar = odefunc_SP1(~,var,miu,m1,m2)
L = 2;
r = L/2;
g = 10;
J = m2*(2*r)^2/12;

th = var(3);
dx = var(4);
dth = var(6);

% slip
d2x = -(2*(J*g*m1*miu*sign(dx) - (r^2*g*m2^2*sin(2*th))/2 - J*r*dth^2*m2*sin(th) - r^3*dth^2*m2^2*sin(th) + J*g*m2*miu*sign(dx) + r^2*g*m1*m2*miu*sign(dx) + r^3*dth^2*m2^2*miu*sign(dx)*cos(th) + r^2*g*m2^2*miu*sign(dx)*cos(th)^2 + J*r*dth^2*m2*miu*sign(dx)*cos(th)))/(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx)*sin(2*th));
d2y = 0;
d2th = -(2*r*m2*(sin(th) - miu*sign(dx)*cos(th))*(r*m2*cos(th)*dth^2 + g*m1 + g*m2))/(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx)*sin(2*th));

dvar(1,1) = var(4);
dvar(2,1) = var(5);
dvar(3,1) = var(6);
dvar(4,1) = d2x;
dvar(5,1) = d2y;
dvar(6,1) = d2th;