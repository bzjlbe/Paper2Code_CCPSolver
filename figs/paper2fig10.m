% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

close all;clear
miu1 = 0.8;                         % miu: Friction coefficient
miu2 = 0.6;
miu3 = 0.4;
rho = 0.3;                          % rho: The parameter of generalized_alpha method 
h = 0.002;                          % h: The time step
num = 300;                          % num: Number of steps
v0 = 2.5;                           % v0: Slider tangential initial velocity
th_slope = pi/12;                   % th_slope: Slope angle of the slope

% CCP_Wang method
sw_on = 0;
    output_Simulation = PM_nonsm_slider(miu1,rho,h,num,v0,th_slope,sw_on);
    time1 = output_Simulation.time;
    F_end1 = output_Simulation.F_end;
    
    output_Simulation = PM_nonsm_slider(miu2,rho,h,num,v0,th_slope,sw_on);
    time2 = output_Simulation.time;
    F_end2 = output_Simulation.F_end;
    
    output_Simulation = PM_nonsm_slider(miu3,rho,h,num,v0,th_slope,sw_on);
    time3 = output_Simulation.time;
    F_end3 = output_Simulation.F_end;

% the proposed method in the article
sw_on = 1;
    output_Simulation = PM_nonsm_slider(miu1,rho,h,num,v0,th_slope,sw_on);
    time4 = output_Simulation.time;
    F_end4 = output_Simulation.F_end;
    
    output_Simulation = PM_nonsm_slider(miu2,rho,h,num,v0,th_slope,sw_on);
    time5 = output_Simulation.time;
    F_end5 = output_Simulation.F_end;
  
    output_Simulation = PM_nonsm_slider(miu3,rho,h,num,v0,th_slope,sw_on);
    time6 = output_Simulation.time;
    F_end6 = output_Simulation.F_end;

    figure(10);
    set(gcf, 'Position', [800 300 1070 420])
    subplot('position',[0.05,0.1,0.42,0.84]);hold on;
        plot(time1,F_end1(3,:),'r',LineWidth=2)
        plot(time2,F_end2(3,:),'-.g',LineWidth=2)
        plot(time3,F_end3(3,:),'--b',LineWidth=2)
        legend('\mu=0.8','\mu=0.6','\mu=0.4'); % ,'location','best' Numerical method,'\mu=0.1'
        xlabel('t/s');ylabel('Normal force/N');ytickformat('%.1f')
        ylim([3,13])
        title('CCP\_Wang');
        % 小图左
        axes_position = [0.28, 0.18, 0.15, 0.4];
        ax_inset1 = axes('position', axes_position);box on
        plot(time1,F_end1(3,:),'r',LineWidth=2);hold on;
        plot(time2,F_end2(3,:),'-.g',LineWidth=2)
        plot(time3,F_end3(3,:),'--b',LineWidth=2)
        ylim([4,10.2])
        xlim([0.0,0.01])
        xtickformat('%.3f')
        ax_inset1.XAxis.Exponent = 0;
        ytickformat('%.1f')
        
    subplot('position',[0.55,0.1,0.42,0.84]);hold on;
        plot(time4,F_end4(3,:),'r',LineWidth=2)
        plot(time5,F_end5(3,:),'-.g',LineWidth=2)
        plot(time6,F_end6(3,:),'--b',LineWidth=2)
        legend('\mu=0.8','\mu=0.6','\mu=0.4'); % ,'location','best' Numerical method,'\mu=0.1'
        xlabel('t/s');ylabel('Normal force/N');ytickformat('%.1f')
        ylim([3,13])
        title('Proposed method');