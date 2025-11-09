% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %


close all;clear
miu1 = 0.1;                                      % miu: Friction coefficient
miu2 = 0.2;
miu3 = 0.3;
rho = 0.3;                                      % rho: The parameter of generalized_alpha method 
h = 0.0001;                                     % h: The time step
num = 4000;                                     % num: Number of steps

% CCP_Wang method
output_Simulation1 = PM_nonsm_Drod(miu1,rho,h,num,0);
time1 = output_Simulation1.time;
q1 = output_Simulation1.q;
F_end1 = output_Simulation1.F_end;

output_Simulation1 = PM_nonsm_Drod(miu2,rho,h,num,0);
time2 = output_Simulation1.time;
q2 = output_Simulation1.q;
F_end2 = output_Simulation1.F_end;

output_Simulation1 = PM_nonsm_Drod(miu3,rho,h,num,0);
time3 = output_Simulation1.time;
q3 = output_Simulation1.q;
F_end3 = output_Simulation1.F_end;

% the proposed method in the article
output_Simulation1 = PM_nonsm_Drod(miu1,rho,h,num,1);
time4 = output_Simulation1.time;
q4 = output_Simulation1.q;
F_end4 = output_Simulation1.F_end;

output_Simulation1 = PM_nonsm_Drod(miu2,rho,h,num,1);
time5 = output_Simulation1.time;
q5 = output_Simulation1.q;
F_end5 = output_Simulation1.F_end;

output_Simulation1 = PM_nonsm_Drod(miu3,rho,h,num,1);
time6 = output_Simulation1.time;
q6 = output_Simulation1.q;
F_end6 = output_Simulation1.F_end;

figure(15);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84])
    plot(time1,F_end1(6,:),'r',LineWidth=2);hold on
    plot(time2,F_end2(6,:),'-.g',LineWidth=2)
    plot(time3,F_end3(6,:),'--b',LineWidth=2)
    xlabel('t/s');ylabel('Normal force F_2/N');
    xtickformat('%.2f')
    ytickformat('%.1f')
    ylim([0,1.8])
    line([0,0.1],[1.4,0.85],Color='black')
    legend('\mu=0.1','\mu=0.2','\mu=0.3'); % ,'location','best' Numerical method
    title('CCP\_Wang');

    % Small picture
    axes_position = [0.08, 0.2, 0.15, 0.3];% left, bottom, width, height
    ax_inset2 = axes('position', axes_position);box on
    % plot(time3,F2_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(6,:),'r',LineWidth=2);hold on
    plot(time2,F_end2(6,:),'-.g',LineWidth=2)
    plot(time3,F_end3(6,:),'--b',LineWidth=2)
    xlim([0.0,0.0005])
    ylim([1.25,1.7])
    ax_inset2.XAxis.Exponent = 0;
    ytickformat('%.1f')
    
subplot('position',[0.55,0.1,0.42,0.84])
    plot(time4,F_end4(6,:),'r',LineWidth=2);hold on
    plot(time5,F_end5(6,:),'-.g',LineWidth=2)
    plot(time6,F_end6(6,:),'--b',LineWidth=2)
    xlabel('t/s');ylabel('Normal force F_2/N');
    xtickformat('%.2f')
    ytickformat('%.1f')
    ylim([0,1.8])
    line([0,0.1],[1.4,0.85],Color='black')
    legend('\mu=0.1','\mu=0.2','\mu=0.3'); % ,'location','best' Numerical method
    title('Proposed method');

    % Small picture
    axes_position = [0.58, 0.2, 0.15, 0.3];% left, bottom, width, height
    ax_inset2 = axes('position', axes_position);box on
    % plot(time3,F2_3,'k',LineWidth=2);hold on
    plot(time4,F_end4(6,:),'r',LineWidth=2);hold on
    plot(time5,F_end5(6,:),'-.g',LineWidth=2)
    plot(time6,F_end6(6,:),'--b',LineWidth=2)
    xlim([0.0,0.0005])
    ylim([1.25,1.7])
    ax_inset2.XAxis.Exponent = 0;
    ytickformat('%.1f')