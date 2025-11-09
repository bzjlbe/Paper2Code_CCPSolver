% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% close all;clear
miu = 0.6;                          % miu: Friction coefficient
rho = 0.3;                          % rho: The parameter of generalized_alpha method 
h = 0.001;                          % h: The time step
num = 600;                          % num: Number of steps
v0 = 2.5;                           % v0: Slider tangential initial velocity
th_slope = pi/12;                   % th_slope: Slope angle of the slope

output_Simulation = PM_nonsm_slider(miu,rho,h,num,v0,th_slope,0);
time = output_Simulation.time;
q = output_Simulation.q;
dq = output_Simulation.dq;

output_Analysis = Analytical_slider(miu,h,num,v0,th_slope);
t_Ana = output_Analysis.time;
diap = output_Analysis.displacement;
velo = output_Analysis.velocity;

figure(3);
% Get the screen size：scnsize = get(0,'ScreenSize')      scnsize = [1,1,1920,1080]
% Set the size and position of the figure displayed on the screen:
% set(gcf, 'Position', [left, bottom, width, height])
set(gcf, 'Position', [400 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84]);
    % Main image, left
    box on;hold on
    plot(t_Ana,diap,'k',LineWidth=2)
    plot(time,q(1,:),'-.r',LineWidth=2)
    
    xlabel('t/s');ylabel('Tangential displacement(m)');ytickformat('%.2f')
    legend('Analytical','CCP\_Wang','location','best','Position',[0.34, 0.7, 0.1, 0.08]);%
       
    rectangle('Position',[0.4, 0.37, 0.15, 0.02], 'EdgeColor','black','Linewidth',2,'LineStyle','--'); % Draw a bounding box in the main image
        % Small graph: Draws the data of the enlarged area in the sub-coordinate system
        axes('position', [0.25 0.3 0.2 0.3]); % [left, bottom, width, height])
        box on;hold on
        plot(t_Ana, diap,'k',LineWidth=2);  
        plot(time,q(1,:),'-.r',LineWidth=2)
        ylim([0.379,0.381])
        xlim([0.40,0.55]);xtickformat('%.2f')
        ytickformat('%.4f')

subplot('position',[0.55,0.1,0.42,0.84]);
    box on;hold on
    plot(t_Ana,velo,'k',LineWidth=2)
    plot(time,dq(1,:),'-.r',LineWidth=2)
    
    xlabel('t/s');ylabel('Tangential velocity(m·s^-1)');ytickformat('%.1f')
    legend('Analytical','CCP\_Wang','location','best');

    rectangle('Position',[0.29, -0.05, 0.04, 0.1], 'EdgeColor','black','Linewidth',2,'LineStyle','--');
        ax_inset = axes('position', [0.75 0.5 0.2 0.3]);
        box on;hold on
        plot(t_Ana,velo,'k',LineWidth=2)
        plot(time,dq(1,:),'-.r',LineWidth=2)
        xlim([0.28,0.35])
        ylim([-0.5e-6,1e-6])
        ax_inset.YAxis.Exponent = -6;
        xtickformat('%.2f')
        ytickformat('%.1f')