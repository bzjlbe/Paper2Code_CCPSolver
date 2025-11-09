% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% close all;
clear;clc
miu = 0.6;                          % miu: Friction coefficient
rho = 0.3;                          % rho: The parameter of generalized_alpha method 
v0 = 2.5;                           % v0: Slider tangential initial velocity
th_slope = pi/12;                   % th_slope: Slope angle of the slope

x = [1e-4,5e-4,1e-3,5e-3,1e-2-1e-3];  % The time step

for i = 1:length(x)
    h = x(i);
    num = round(0.6/h);
    output_Analysis = Analytical_slider(miu,h,num,v0,th_slope);
    output_Simulation1 = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,0);
    output_Simulation2 = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,1);
    err_disp_Wang(i) = sum(abs(output_Analysis.displacement - output_Simulation1.q(1,:)))/sum(abs(output_Analysis.displacement));
    err_velo_Wang(i) = sum(abs(output_Analysis.velocity - output_Simulation1.dq(1,:)))/sum(abs(output_Analysis.velocity));
    err_disp_proposed(i) = sum(abs(output_Analysis.displacement - output_Simulation2.q(1,:)))/sum(abs(output_Analysis.displacement));
    err_velo_proposed(i) = sum(abs(output_Analysis.velocity - output_Simulation2.dq(1,:)))/sum(abs(output_Analysis.velocity));
end
figure(12)
set(gcf, 'Position', [400 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84]);
    loglog(x,err_disp_Wang,'-or',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_disp_proposed,'-xg',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential displacement error')
    legend('disp\_Wang','disp\_proposed')% ,Location='best'


subplot('position',[0.55,0.1,0.42,0.84]);
    loglog(x,err_velo_Wang,'--pb',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_velo_proposed,'--sk',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential velocity error')
    legend('velo\_Wang','velo\_proposed',Location='best')