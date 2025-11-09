% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% close all;
clear;clc
miu = 0.1;                                      % miu: Friction coefficient
rho = 0.3;                                      % rho: The parameter of generalized_alpha method 
h = 0.0001;                                     % h: The time step
num = 4000;                                     % num: Number of steps

x = [1e-5,1e-4,1e-3,1e-2];
y = 100*[1,10,100,1000];

% Reference Solution
output_Analysis = ODE_Drod(miu,1e-7);
for i = 1:length(x)
    h = x(i);
    num = round(0.4/h);
    
    output_Simulation1 = PM_nonsm_Drod(miu,rho,h,num,0);        % CCP_Wang method
    output_Simulation2 = PM_nonsm_Drod(miu,rho,h,num,1);        % the proposed method in the article

    err_disp_Wang(i) = sum(abs(output_Analysis.th(1:y(i):end) - output_Simulation1.th))/sum(abs(output_Analysis.th(1:y(i):end)));
    err_velo_Wang(i) = sum(abs(output_Analysis.dth(1:y(i):end) - output_Simulation1.dth))/sum(abs(output_Analysis.dth(1:y(i):end)));
    err_disp_proposed(i) = sum(abs(output_Analysis.th(1:y(i):end) - output_Simulation2.th))/sum(abs(output_Analysis.th(1:y(i):end)));
    err_velo_proposed(i) = sum(abs(output_Analysis.dth(1:y(i):end) - output_Simulation2.dth))/sum(abs(output_Analysis.dth(1:y(i):end)));
end

figure(15)
set(gcf, 'Position', [400 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84]);
    loglog(x,err_disp_Wang,'-or',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_disp_proposed,'-xg',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential displacement error of the right end')
    legend('disp\_Wang','disp\_proposed')% ,Location='best'

subplot('position',[0.55,0.1,0.42,0.84]);
    loglog(x,err_velo_Wang,'--pb',LineWidth=2,MarkerSize=10);hold on
    loglog(x,err_velo_proposed,'--sk',LineWidth=2,MarkerSize=10)
    xlabel('time step size(s)')
    ylabel('Tangential velocity error of the right end')
    legend('velo\_Wang','velo\_proposed',Location='best')