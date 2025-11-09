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
    output_Analysis = Analytical_slider(miu1,h,num,v0,th_slope);
    output_Simulation = PM_nonsm_slider(miu1,rho,h,num,v0,th_slope,sw_on);
    time4 = output_Simulation.time;
    F_end4 = output_Simulation.F_end;

    output_Simulation = PM_nonsm_slider(miu2,rho,h,num,v0,th_slope,sw_on);
    time5 = output_Simulation.time;
    F_end5 = output_Simulation.F_end;

    output_Simulation = PM_nonsm_slider(miu3,rho,h,num,v0,th_slope,sw_on);
    time6 = output_Simulation.time;
    F_end6 = output_Simulation.F_end;

    figure(11);
    set(gcf, 'Position', [800 300 1070 420])
    subplot('position',[0.05,0.1,0.42,0.84])
        plot(time1,-F_end1(1,:),'r',LineWidth=2);hold on;
        plot(time2,-F_end2(1,:),'-.g',LineWidth=2)
        plot(time3,-F_end3(1,:),'--b',LineWidth=2)
        legend('\mu=0.8','\mu=0.6','\mu=0.4');  % ,'location','best'
        xlabel('t/s');ylabel('Friction/N');
        ylim([-8.5,9])
        title('CCP\_Wang');
        
    subplot('position',[0.55,0.1,0.42,0.84])
        plot(time4,-F_end4(1,:),'r',LineWidth=2);hold on;
        plot(time5,-F_end5(1,:),'-.g',LineWidth=2)
        plot(time6,-F_end6(1,:),'--b',LineWidth=2)
        legend('\mu=0.8','\mu=0.6','\mu=0.4');  % ,'location','best'
        xlabel('t/s');ylabel('Friction/N');
        ylim([-8.5,9])
        title('Proposed method');