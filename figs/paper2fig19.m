% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %


% close all;
clear;clc
miu = 0.1;                          % Friction coefficient
th0 = pi/2;                         % Initial value of rod angle

dx0 = -1;                           % Initial horizontal velocity of the slider
m1 = 3;                             % Slider quality
m2 = 1;                             % Rod quality

L = 2;                              % the length of the rod
g = 10;                             % Gravity
J = m2*L^2/12;                      % Inertia

rho = 0.3;                          % The parameter of generalized_alpha method 
h = 0.001;                          % The time step
num = 4000;                         % Number of steps

% Reference Solution
output_ana = ODE_SP(miu,th0,dx0,m1,m2,h);
time1 = output_ana.time;
th1 = output_ana.th;
dth1 = output_ana.dth;
d2th1 = output_ana.d2th;
x1 = output_ana.x;
dx1 = output_ana.dx;
d2x1 = output_ana.d2x;
Fx1 = output_ana.Fx;
Fz1 = output_ana.Fz;

Func_x1 = m1*d2x1+m2*(d2x1+L*d2th1.*cos(th1)-L*dth1.^2.*sin(th1)) - Fx1;
Func_y1 = m2*(d2th1*L.*sin(th1)+dth1.^2*L.*cos(th1)) - Fz1 + (m1+m2)*g;
Func_th1 = J*d2th1 + L*m2*(L*d2th1 + d2x1.*cos(th1) + g*sin(th1));

E1 = 0.5*m1*dx1.^2 + 0.5*m2*((dx1+L*dth1.*cos(th1)).^2+(L*dth1.*sin(th1)).^2) + 0.5*J*(L*dth1).^2 - m2*g*L*cos(th1);
E1mgh =  - m2*g*L*cos(th1);
E1v =  0.5*m1*dx1.^2 + 0.5*m2*((dx1+L*dth1.*cos(th1)).^2+(L*dth1.*sin(th1)).^2) + 0.5*J*(L*dth1).^2;

% CCP_Wang method
output_Wang = PM_nonsm_SP(miu,th0,dx0,m1,m2,rho,h,num,0);
time2 = output_Wang.time;
th2 = output_Wang.th;
dth2 = output_Wang.dth;
d2th2 = output_Wang.d2th;
x2 = output_Wang.x;
dx2 = output_Wang.dx;
d2x2 = output_Wang.d2x;
Fx2 = output_Wang.Fx;
Fz2 = output_Wang.Fz;

Func_x2 = m1*d2x2+m2*(d2x2+L*d2th2.*cos(th2)-L*dth2.^2.*sin(th2)) - Fx2;
Func_y2 = m2*(d2th2*L.*sin(th2)+dth2.^2*L.*cos(th2)) - Fz2 + (m1+m2)*g;
Func_th2 = J*d2th2 + L*m2*(L*d2th2 + d2x2.*cos(th2) + g*sin(th2));

E2 = 0.5*m1*dx2.^2 + 0.5*m2*((dx2+L*dth2.*cos(th2)).^2+(L*dth2.*sin(th2)).^2) + 0.5*J*(L*dth2).^2 - m2*g*L*cos(th2);
E2mgh =  - m2*g*L*cos(th2);
E2v = 0.5*m1*dx2.^2 + 0.5*m2*((dx2+L*dth2.*cos(th2)).^2+(L*dth2.*sin(th2)).^2) + 0.5*J*(L*dth2).^2;

% the proposed method in the article
output = PM_nonsm_SP(miu,th0,dx0,m1,m2,rho,h,num,1);
time3 = output.time;
th3 = output.th;
dth3 = output.dth;
d2th3 = output.d2th;
x3 = output.x;
dx3 = output.dx;
d2x3 = output.d2x;
Fx3 = output.Fx;
Fz3 = output.Fz;

E3 = 0.5*m1*dx3.^2 + 0.5*m2*((dx3+L*dth3.*cos(th3)).^2+(L*dth3.*sin(th3)).^2) + 0.5*J*(L*dth3).^2 - m2*g*L*cos(th3);
E3mgh = - m2*g*L*cos(th3);
E3v = 0.5*m1*dx3.^2 + 0.5*m2*((dx3+L*dth3.*cos(th3)).^2+(L*dth3.*sin(th3)).^2) + 0.5*J*(L*dth3).^2;

figure(18);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.06,0.1,0.42,0.84])
    plot(time1,x1,'k',LineWidth=2);hold on
    plot(time2,x2,'-.r',LineWidth=2)
    plot(time3,x3,'--g',LineWidth=2)
    xlabel('t/s');ylabel('x/m');
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method,'CCP\_Wang'

subplot('position',[0.56,0.1,0.42,0.84])
    plot(time1,th1,'k',LineWidth=2);hold on
    plot(time2,th2,'-.r',LineWidth=2)
    plot(time3,th3,'--g',LineWidth=2)
    xlabel('t/s');ylabel('\theta/rad');
    ylim([-1.8,2])
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method'CCP\_Wang',
