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
x = [1e-4,5e-4,1e-3,5e-3,1e-2];     % The time step
y = 1000*[1,5,10,50,100];

output_Analysis = ODE_SP(miu,th0,dx0,m1,m2,1e-7);
for i = 1:length(x)
    h = x(i);
    num = round(4/h);
    output_Simulation1 = PM_nonsm_SP(miu,th0,dx0,m1,m2,rho,h,num,0);
    output_Simulation2 = PM_nonsm_SP(miu,th0,dx0,m1,m2,rho,h,num,1);

    err_disp_Wang(i) = sum(abs(output_Analysis.x(1:y(i):end) - output_Simulation1.x.'))/sum(abs(output_Analysis.x(1:y(i):end)));
    err_velo_Wang(i) = sum(abs(output_Analysis.dx(1:y(i):end) - output_Simulation1.dx.'))/sum(abs(output_Analysis.dx(1:y(i):end)));
    err_disp_proposed(i) = sum(abs(output_Analysis.x(1:y(i):end) - output_Simulation2.x.'))/sum(abs(output_Analysis.x(1:y(i):end)));
    err_velo_proposed(i) = sum(abs(output_Analysis.dx(1:y(i):end) - output_Simulation2.dx.'))/sum(abs(output_Analysis.dx(1:y(i):end)));
end

figure(23)
set(gcf, 'Position', [400 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84]);
    loglog(x,err_disp_Wang,'-or',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_disp_proposed,'-xg',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential displacement error of the slider')
    legend('disp\_Wang','disp\_proposed')% ,Location='best'

subplot('position',[0.55,0.1,0.42,0.84]);
    loglog(x,err_velo_Wang,'--pb',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_velo_proposed,'--sk',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential velocity error of the slider')
    legend('velo\_Wang','velo\_proposed',Location='best')