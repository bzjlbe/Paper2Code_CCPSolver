function output = ODE_Drod(miu,h)
L = 1;
r = 0.002;
m = 0.402;
J = m*L*L/12;
g = 9.81;

refine = 4;
options = odeset('Events',@(t,y) events(t,y,miu),'Refine',refine,'RelTol',1e-6); % 'OutputFcn',@odeplot,'OutputSel',1,

tstart = 0;
tfinal = 0.4;

% Simulation Section 1: Contact on both sides
th0 = pi/4;
dth0 = 0;
var_th0 = [th0;dth0];
tic
[t1,var_th] = ode45(@(t,var) odefunc1(t,var,miu),[tstart:h:tfinal],var_th0,options);

th = var_th(:,1);
dth = var_th(:,2);
x1 = L*cos(th)/2 + r*sin(th);
dx1 = - L*sin(th).*dth/2 + r*cos(th).*dth;
y1 = L*sin(th)/2 + r*cos(th);
dy1 = L*cos(th).*dth/2 - r*sin(th).*dth;

A = (- m*L^2*miu^2 + m*L^2 - 4*m*cos(2*th)*L*miu*r - 4*m*sin(2*th)*L*r + 4*m*miu^2*r^2 + 4*J*miu^2 + 4*m*r^2 + 4*J)/(4*(miu^2 + 1));
B = (L*m*(L*miu - 2*r*cos(2*th) + 2*miu*r*sin(2*th)))/(2*(miu^2 + 1));
C = (g*m*(2*r*sin(th) - L*cos(th) + 2*L*miu*sin(th) + L*miu^2*cos(th) + 2*miu^2*r*sin(th)))/(2*(miu^2 + 1));
d2th1 = (C-B.*dth.^2)./A;

F1_1 = (m*(g - d2th1*r.*sin(th) - (L*dth.^2.*sin(th))/2 - dth.^2*r.*cos(th) + (L*d2th1.*cos(th))/2 + (L*d2th1*miu.*sin(th))/2 - d2th1*miu*r.*cos(th) + (L*dth.^2*miu.*cos(th))/2 + dth.^2*miu*r.*sin(th)))/(miu^2 + 1);
F2_1 = -(m*((L*dth.^2.*cos(th))/2 - d2th1*r.*cos(th) - g*miu + dth.^2*r.*sin(th) + (L*d2th1.*sin(th))/2 - (L*d2th1*miu.*cos(th))/2 + d2th1*miu*r.*sin(th) + (L*dth.^2*miu.*sin(th))/2 + dth.^2*miu*r.*cos(th)))/(miu^2 + 1);

% Simulation Section 2: Left end detached
nt = floor(t1(end)/h)+1;
options = odeset(options,'InitialStep',t1(nt)-t1(nt-refine),'MaxStep',t1(nt)-t1(1),'RelTol',1e-6);

tstart = t1(nt);

th0 = var_th(nt,1);
dth0 = var_th(nt,2);

x0 = L*cos(th0)/2 + r*sin(th0);
dx0 = - L*sin(th0)*dth0/2 + r*cos(th0)*dth0;

var_th_x0(1:2,1) = [th0;x0];
var_th_x0(3:4,1) = [dth0;dx0];


if tstart ~= tfinal
    [t2,var_th_x] = ode45(@(t,y) odefunc2(t,y,miu),tstart:h:tfinal,var_th_x0,options);
    toc
    th = var_th_x(:,1);
    dth = var_th_x(:,3);
    y2 = L*sin(th)/2 + r*cos(th);
    dy2 = L*cos(th).*dth/2 - r*sin(th).*dth;
    
    A = J + (L*m*(L*cos(th) - 2*r*sin(th)).*(cos(th) - miu*sin(th)))/4 - (m*r*(L*cos(th) - 2*r*sin(th)).*(sin(th) + miu*cos(th)))/2;
    B = (m*(L*sin(th) + 2*r*cos(th)).*(2*r*sin(th) - L*cos(th) + L*miu*sin(th) + 2*miu*r*cos(th)))/4;
    C = (g*m*(2*r*sin(th) - L*cos(th) + L*miu*sin(th) + 2*miu*r*cos(th)))/2;
    
    d2th2 = (C-B.*dth.^2)./A;
    
    F1_2 = -m*(d2th2*r.*sin(th) - g + (L*dth.^2.*sin(th))/2 + dth.^2*r.*cos(th) - (L*d2th2.*cos(th))/2);
    F2_2 = 0.*F1_2;

    output.time = [t1(1:nt-1);t2].';
    output.th = [var_th(1:nt-1,1);var_th_x(:,1)].';
    output.dth = [var_th(1:nt-1,2);var_th_x(:,3)].';
    output.d2th = [d2th1(1:nt-1);d2th2].';
    output.x = [x1(1:nt-1);var_th_x(:,2)].';
    output.dx = [dx1(1:nt-1);var_th_x(:,4)].';
    output.y = [y1(1:nt-1);y2].';
    output.dy = [dy1(1:nt-1);dy2].';
    output.F1 = [F1_1(1:nt-1);F1_2].';
    output.F2 = [F2_1(1:nt-1);F2_2].';

else
output.time = t1.';
output.th = var_th(:,1).';
output.dth = var_th(:,2).';
output.d2th = d2th1.';
output.x = x1.';
output.dx = dx1.';
output.y = y1.';
output.dy = dy1.';
output.F1 = F1_1.';
output.F2 = F2_1.';
end


% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y,miu)
L = 1;
r = 0.002;
m = 0.402;
g = 9.81;

dy = odefunc1(t,y,miu);

th = y(1);
dth = y(2);
d2th = dy(2);

value = -(m*((L*dth^2*cos(th))/2 - d2th*r*cos(th) - g*miu + dth^2*r*sin(th) + (L*d2th*sin(th))/2 - (L*d2th*miu*cos(th))/2 + d2th*miu*r*sin(th) + (L*dth^2*miu*sin(th))/2 + dth^2*miu*r*cos(th)))/(miu^2 + 1);
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
