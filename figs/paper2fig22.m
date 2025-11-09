% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% close all;
clear;clc
miu = 0.8;                          % Friction coefficient
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

% the proposed method in the article
output = PM_nonsm_SP(miu,th0,dx0,m1,m2,rho,h,num,1);
time3 = output.time;
th3 = output.th;
dth3 = output_Wang.dth;
d2th3 = output_Wang.d2th;
x3 = output.x;
dx3 = output.dx;
d2x3 = output.d2x;
Fx3 = output.Fx;
Fz3 = output.Fz;

E3 = 0.5*m1*dx3.^2 + 0.5*m2*((dx3+L*dth3.*cos(th3)).^2+(L*dth3.*sin(th3)).^2) + 0.5*J*(L*dth3).^2 - m2*g*L*cos(th3);

figure(22);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.06,0.1,0.42,0.84])
    plot(time1,Fz1,'k',LineWidth=2);hold on
    plot(time2,Fz2,'-.r',LineWidth=2)
    plot(time3,Fz3,'--g',LineWidth=2)
    xlabel('t/s');ylabel('Normal force/N');
    ylim([20,60])
    line([0,1],[32.5,28],Color='black')
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method

    axes_position = [0.1, 0.14, 0.14, 0.2]; % left, bottom, width, height
    ax_inset1 = axes('position', axes_position);box on
    plot(time1,Fz1,'k',LineWidth=2);hold on
    plot(time2,Fz2,'-.r',LineWidth=2)
    plot(time3,Fz3,'--g',LineWidth=2)
    xlim([0.0,0.05])
    ylim([19,37])
    ax_inset1.XAxis.Exponent = 0;
    xtickformat('%.2f')

subplot('position',[0.56,0.1,0.42,0.84])
    plot(time1,Fx1,'k',LineWidth=2);hold on
    plot(time2,Fx2,'-.r',LineWidth=2)
    plot(time3,Fx3,'--g',LineWidth=2)
    xlabel('t/s');ylabel('Friction/N');
    line([0.2,1.6],[25,21],Color='black')
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method

    axes_position = [0.66, 0.7, 0.14, 0.2];
    ax_inset1 = axes('position', axes_position);box on 
    plot(time1,Fx1,'k',LineWidth=2);hold on
    plot(time2,Fx2,'-.r',LineWidth=2)
    plot(time3,Fx3,'--g',LineWidth=2)
    xlim([0.0,0.05])
    ylim([15,30])
    ax_inset1.XAxis.Exponent = 0;
    xtickformat('%.2f')