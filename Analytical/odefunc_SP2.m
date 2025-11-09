% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

function dvar = odefunc_SP2(~,var,m2)
L = 2;
r = L/2;
g = 10;
J = m2*(2*r)^2/12;

th = var(3);

% stick
d2x = 0;
d2y = 0;
d2th = -(r*g*m2*sin(th))/(m2*r^2 + J);

dvar(1,1) = var(4);
dvar(2,1) = var(5);
dvar(3,1) = var(6);
dvar(4,1) = d2x;
dvar(5,1) = d2y;
dvar(6,1) = d2th;