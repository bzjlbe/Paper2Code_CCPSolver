% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

close all;clear
miu = 0.6;                          % miu: Friction coefficient
rho = 0.3;                          % rho: The parameter of generalized_alpha method 
h = 0.001;                          % h: The time step
num = 600;                          % num: Number of steps
v0 = 2.5;                           % v0: Slider tangential initial velocity
th_slope = pi/12;                   % th_slope: Slope angle of the slope

output_Simulation = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,1);
time = output_Simulation.time;
dmvt = output_Simulation.dmvt;

output_Simulation = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,0);
dq = output_Simulation.dq;

figure(6);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84])
    plot(time,dq(2,:),'k',LineWidth=2);hold on
    plot(time,dmvt,'--k',LineWidth=2)
    ylim([-0.0035,0.0035])
    title('Velocity errors')
    legend('Normal velocity error','Predictive tangential velocity error'); % ,'location','best' Numerical method
    xlabel('t/s');ylabel('Velocity error/(m/s)');ytickformat('%.1f')

subplot('position',[0.55,0.1,0.42,0.84])
    plot(time,dq(2,:) + dmvt,'k',LineWidth=2)
    title('The sum of the two errors');
    xlabel('t/s');ylabel('Sum of velocity errors/N');
    % ylim([-4,4])
