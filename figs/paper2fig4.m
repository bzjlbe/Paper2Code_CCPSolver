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

output_Simulation = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,0);
time = output_Simulation.time;
dq = output_Simulation.dq;
F_end = output_Simulation.F_end;

output_Analysis = Analytical_slider(miu,h,num,v0,th_slope);
t_Ana = output_Analysis.time;
diap = output_Analysis.displacement;
velo = output_Analysis.velocity;
F_Ana = output_Analysis.F_Analytical;

figure(4);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84])
    plot(t_Ana,F_Ana(3,:),'k',LineWidth=2);hold on;
    plot(time,F_end(3,:),'-.r',LineWidth=2)
    line([0.04,0.01],[6,7],Color='black')
    line([0.34,0.31],[6,9],Color='black')
    legend('Analytical','CCP\_Wang'); % ,'location','best' Numerical method
    xlabel('t/s');ylabel('Normal force/N');
    % ytickformat('%.1f')
    ylim([3,13])
    % title('Normal force');
    % Small picture, left
    axes_position = [0.09, 0.18, 0.15, 0.4]; % left, bottom, width, height
    ax_inset1 = axes('position', axes_position);box on
    plot(t_Ana,F_Ana(3,:),'k',LineWidth=2);hold on;
    plot(time,F_end(3,:),'-.r',LineWidth=2)
    ylim([5.5,10.5])
    xlim([0.0,0.006])
    ax_inset1.XAxis.Exponent = 0;
    % ytickformat('%.1f')
    % Small picture, right
    axes_position = [0.3, 0.18, 0.15, 0.4];
    axes('position', axes_position);box on 
    plot(t_Ana,F_Ana(3,:),'k',LineWidth=2);hold on; % 解析(1:25:end)
    plot(time,F_end(3,:),'-.r',LineWidth=2)
    ylim([8.5,13])
    xlim([0.303,0.310])
    xtickformat('%.3f')
    % ytickformat('%.1f')
    
subplot('position',[0.55,0.1,0.42,0.84])
    plot(t_Ana,-F_Ana(1,:),'k',LineWidth=2);hold on 
    plot(time,-F_end(1,:),'-.r',LineWidth=2)
    line([0.15,0.01],[3.8,5],Color='black')
    line([0.25,0.3],[-3.8,-2.4],Color='black')
    legend('Analytical','CCP\_Wang');  % ,'location','best'
    xlabel('t/s');ylabel('Friction/N');
    ylim([-5.5,6.5])
    % title('Friction');
    % Small picture, top
    axes_position = [0.58, 0.5, 0.15, 0.25]; 
    ax_inset4 = axes('position', axes_position);box on 
    plot(t_Ana,-F_Ana(1,:),'k',LineWidth=2);hold on  
    plot(time,-F_end(1,:),'-.r',LineWidth=2)
    ylim([3.0,6.1])
    xlim([0.0,0.005])
    ax_inset4.XAxis.Exponent = 0;
    % ytickformat('%.1f')
    % Small picture, below
    axes_position = [0.58, 0.17, 0.15, 0.25];
    axes('position', axes_position);box on 
    plot(t_Ana,-F_Ana(1,:),'k',LineWidth=2);hold on 
    plot(time,-F_end(1,:),'-.r',LineWidth=2)
    ylim([-5.5,-2.0])
    xlim([0.303,0.308])
    xtickformat('%.3f')
    % ytickformat('%.1f')