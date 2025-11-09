% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% close all;clear
miu = 0.1;                                      % miu: Friction coefficient
rho = 0.3;                                      % rho: The parameter of generalized_alpha method 
h = 0.0001;                                     % h: The time step
num = 4000;                                     % num: Number of steps

% CCP_Wang method
output_Simulation1 = PM_nonsm_Drod(miu,rho,h,num,0);
time1 = output_Simulation1.time;
q1 = output_Simulation1.q;
F_end1 = output_Simulation1.F_end;

% the proposed method in the article
output_Simulation2 = PM_nonsm_Drod(miu,rho,h,num,1);
time2 = output_Simulation2.time;
q2 = output_Simulation2.q;
F_end2 = output_Simulation2.F_end;

% Reference Solution
output_Analysis = ODE_Drod(miu,h);
time3 = output_Analysis.time;
th3 = output_Analysis.th;
F1_3 = output_Analysis.F1;
F2_3 = output_Analysis.F2;

figure(14);
set(gcf, 'Position', [800 300 1070 420])
subplot('position',[0.05,0.1,0.42,0.84])
    plot(time3,F2_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(6,:),'-.r',LineWidth=2)
    plot(time2,F_end2(6,:),'--g',LineWidth=2)
    xlabel('t/s');ylabel('Normal force F_2/N');
    xtickformat('%.2f')
    ytickformat('%.1f')
    ylim([0,1.8])
    line([0,0.1],[1.4,0.85],Color='black')
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method

    % Small picture
    axes_position = [0.08, 0.2, 0.15, 0.3];% left, bottom, width, height
    ax_inset2 = axes('position', axes_position);box on
    plot(time3,F2_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(6,:),'-.r',LineWidth=2)
    plot(time2,F_end2(6,:),'--g',LineWidth=2)
    xlim([0.0,0.001])
    ylim([1.3,1.7])
    ax_inset2.XAxis.Exponent = 0;
    ytickformat('%.1f')
    
subplot('position',[0.55,0.1,0.42,0.84])
    plot(time3,F1_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(3,:),'-.r',LineWidth=2)
    plot(time2,F_end2(3,:),'--g',LineWidth=2)
    xlabel('t/s');ylabel('Normal force F_1/N');
    xtickformat('%.2f')
    ytickformat('%.1f')
    ylim([0.5,3])
    line([0,0.1],[2.6,1.7],Color='black')
    line([0.268,0.325],[1.0,1.4],Color='black')
    legend('Reference','CCP\_Wang','Proposed method'); % ,'location','best' Numerical method

    % Small picture, left
    axes_position = [0.58, 0.2, 0.15, 0.3];% left, bottom, width, height
    ax_inset1 = axes('position', axes_position);box on
    plot(time3,F1_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(3,:),'-.r',LineWidth=2)
    plot(time2,F_end2(3,:),'--g',LineWidth=2)
    xlim([0.0,0.001])
    ylim([2.4,3.0])
    ax_inset1.XAxis.Exponent = 0;
    ytickformat('%.1f')
    % Small picture, right
    axes_position = [0.81, 0.46, 0.15, 0.3];% left, bottom, width, height
    axes('position', axes_position);box on
    plot(time3,F1_3,'k',LineWidth=2);hold on
    plot(time1,F_end1(3,:),'-.r',LineWidth=2)
    plot(time2,F_end2(3,:),'--g',LineWidth=2)
    xlim([0.268,0.2705])
    ylim([1.005,1.015])
    ytickformat('%.3f')
    xtickformat('%.4f')