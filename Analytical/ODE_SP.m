% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

function output = ODE_SP(miu,th0,dx0,m1,m2,h)
L = 2;
r = L/2;
g = 10;
J = m2*(2*r)^2/12;

refine = 4;
options = odeset('Events',@(t,y) events(t,y,miu,m1,m2),'Refine',refine,'RelTol',1e-6);

tstart = 0;
tfinal = 4;

% Simulation Section 1: Slip
x0 = 0;y0 = 0;dy0 = 0;dth0 = 0;
var0 = [x0;y0;th0;dx0;dy0;dth0];
[t1,var1] = ode45(@(t,var) odefunc_SP1(t,var,miu,m1,m2),[tstart:h:tfinal],var0,options);

th1 = var1(:,3);
dx1 = var1(:,4);
dth1 = var1(:,6);

d2x1 = -(2*(J*g*m1*miu*sign(dx1) - (r^2*g*m2^2*sin(2*th1))/2 - J*r*dth1.^2*m2.*sin(th1) - r^3*dth1.^2*m2^2.*sin(th1) + J*g*m2*miu*sign(dx1) + r^2*g*m1*m2*miu*sign(dx1) + r^3*dth1.^2*m2^2*miu.*sign(dx1).*cos(th1) + r^2*g*m2^2*miu*sign(dx1).*cos(th1).^2 + J*r*dth1.^2*m2*miu.*sign(dx1).*cos(th1)))./(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th1) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx1).*sin(2*th1));
d2y1 = 0.*t1;
d2th1 = -(2*r*m2*(sin(th1) - miu*sign(dx1).*cos(th1)).*(r*m2*cos(th1).*dth1.^2 + g*m1 + g*m2))./(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th1) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx1).*sin(2*th1));
Fx1 = -(2*miu*sign(dx1)*(m1*m2*r^2 + J*m1 + J*m2).*(r*m2*cos(th1).*dth1.^2 + g*m1 + g*m2))./(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th1) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx1).*sin(2*th1));
Fz1 = (2*(m1*m2*r^2 + J*m1 + J*m2)*(r*m2*cos(th1).*dth1.^2 + g*m1 + g*m2))./(2*J*m1 + 2*J*m2 + r^2*m2^2 - r^2*m2^2*cos(2*th1) + 2*r^2*m1*m2 - r^2*m2^2*miu*sign(dx1).*sin(2*th1));

if norm(t1(end) - 4) > 1e-7
    % Simulation Section 2: Stick
    nt = floor(t1(end)/h)+1;
    options = odeset(options,'InitialStep',t1(nt)-t1(nt-refine),'MaxStep',t1(nt)-t1(1),'RelTol',1e-6);
    tstart = t1(nt);
    
    var0 = var1(end,:).';
    
    [t2,var2] = ode45(@(t,var) odefunc_SP2(t,var,miu,m1,m2),[tstart:h:tfinal],var0,options);
    
    th2 = var2(:,3);
    dth2 = var2(:,6);
    
    d2x2 = 0.*t2;
    d2y2 = 0.*t2;
    d2th2 = -(r*g*m2*sin(th2))/(m2*r^2 + J);
    Fx2 = -(r*m2*(2*m2*sin(th2)*r^2.*dth2.^2 + g*m2*sin(2*th2)*r + 2*J*sin(th2).*dth2.^2))/(2*(m2*r^2 + J));
    Fz2 = (r^3*dth2.^2*m2^2.*cos(th2) + g*r^2*m2^2*cos(th2).^2 + g*m1*r^2*m2 + J*r*dth2.^2*m2.*cos(th2) + J*g*m2 + J*g*m1)/(m2*r^2 + J);

    time = [t1(1:nt);t2];
    x = [var1(1:nt,1);var2(:,1)];
    y = [var1(1:nt,2);var2(:,2)];
    th2 = [var1(1:nt,3);var2(:,3)];
    dx2 = [var1(1:nt,4);var2(:,4)];
    dy = [var1(1:nt,5);var2(:,5)];
    dth2 = [var1(1:nt,6);var2(:,6)];
    
    d2x = [d2x1(1:nt);d2x2];
    d2y = [d2y1(1:nt);d2y2];
    d2th = [d2th1(1:nt);d2th2];
    Fx = [Fx1(1:nt);Fx2];
    Fz = [Fz1(1:nt);Fz2];
else
    time = t1;
    x = var1(:,1);
    y = var1(:,2);
    th2 = var1(:,3);
    dx2 = var1(:,4);
    dy = var1(:,5);
    dth2 = var1(:,6);
    
    d2x = d2x1;
    d2y = d2y1;
    d2th = d2th1;
    Fx = Fx1;
    Fz = Fz1;
end

output.time = time;
output.th = th2;
output.x = x;
output.y = y;
output.dth = dth2;
output.dx = dx2;
output.dy = dy;
output.d2th = d2th;
output.d2x = d2x;
output.d2y = d2y;
output.Fx = Fx;
output.Fz = Fz;


function [value,isterminal,direction] = events(~,var)
dx = var(4);
value =  abs(dx)-1e-6;
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
