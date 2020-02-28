
function plot_periodic_results(case_sel, Results)
        
    path = sprintf('figures\\case %d\\', case_sel);
    
    fig1=figure(1);
    set(gca,'FontSize',11)    
    plot(Results(1,:,1),'ro','LineWidth',1)
    hold on;
    plot(Results(2,:,1),'o','LineWidth',1,'color',[0,0,0]) 
    plot(Results(3,:,1),'o','LineWidth',1,'color',[0,0,0.6])
    plot(Results(4,:,1),'o','LineWidth',1,'color',[0.5,0,0.6])
    plot(Results(5,:,1),'o','LineWidth',1,'color',[0.6,0,0])
    plot(Results(7,:,1),'o','LineWidth',1,'color',[0.6,0.6,0])
    plot(Results(10,:,1),'o','LineWidth',1,'color',[0,0.6,0])
    plot(Results(14,:,1),'o','LineWidth',1,'color',[0,0.5,0.8])    
    xlabel('Active Days','fontsize',12)
    ylabel('Total Drug (mg/m^2)','fontsize',12)
    % ylim([0 35])
    legend('Continuous','2-day period','3-day period','4-day period'...
        ,'5-day period','7-day period','10-day period','14-day period');
    title('Comparison of drug dosage', 'fontsize', 12)
    
    saveas(fig1, strcat(path,'comp_v_total'), 'fig');
    print(fig1,'-dpng', strcat(path,'comp_v_total.png'), '-r150');


    fig2=figure(2);
    set(gca,'FontSize',11)
    plot(Results(1,:,2),'ro','LineWidth',1)
    hold on;
    plot(Results(2,:,2),'o','LineWidth',1,'color',[0,0,0]) 
    plot(Results(3,:,2),'o','LineWidth',1,'color',[0,0,0.6])
    plot(Results(4,:,2),'o','LineWidth',1,'color',[0.5,0,0.6])
    plot(Results(5,:,2),'o','LineWidth',1,'color',[0.6,0,0])
    plot(Results(7,:,2),'o','LineWidth',1,'color',[0.6,0.6,0])
    plot(Results(10,:,2),'o','LineWidth',1,'color',[0,0.6,0])
    plot(Results(14,:,2),'o','LineWidth',1,'color',[0,0.5,0.8])
    xlabel('Active Days','fontsize',12)
    ylabel('Maximum drug concentration (mg/L)','fontsize',12)
    % ylim([0 4.5]) % ylim([0.55 1.1])
    legend('Continuous','2-day period','3-day period','4-day period'...
        ,'5-day period','7-day period','10-day period','14-day period');
    title('Comparison of drug concentration', 'fontsize', 12)
    
    saveas(fig2, strcat(path,'comp_u_max'), 'fig');
    print(fig2,'-dpng', strcat(path,'comp_u_max.png'));


    fig3=figure(3);
    set(gca,'FontSize',11)
    plot(Results(1,:,3),'ro','LineWidth',1)
    hold on;
    plot(Results(2,:,3),'o','LineWidth',1,'color',[0,0,0]) 
    plot(Results(3,:,3),'o','LineWidth',1,'color',[0,0,0.6])
    plot(Results(4,:,3),'o','LineWidth',1,'color',[0.5,0,0.6])
    plot(Results(5,:,3),'o','LineWidth',1,'color',[0.6,0,0])
    plot(Results(7,:,3),'o','LineWidth',1,'color',[0.6,0.6,0])
    plot(Results(10,:,3),'o','LineWidth',1,'color',[0,0.6,0])
    plot(Results(14,:,3),'o','LineWidth',1,'color',[0,0.5,0.8])
    xlabel('Active Days','fontsize',12)
    ylabel('Time (Days)','fontsize',12)
%     ylim([0 95]) % 
%     ylim([0 180])
    legend('Continuous','2-day period','3-day period','4-day period'...
        ,'5-day period','7-day period','10-day period','14-day period');
    title('Comparison of tumor eradication time', 'fontsize', 12)
    
    saveas(fig3, strcat(path,'comp_t_zero'), 'fig');
    print(fig3,'-dpng', strcat(path,'comp_t_zero.png'), '-r150');


    fig4=figure(4);
    set(gca,'FontSize',11)
    plot(Results(1,:,4),'ro','LineWidth',1)
    hold on;
    plot(Results(2,:,4),'o','LineWidth',1,'color',[0,0,0]) 
    plot(Results(3,:,4),'o','LineWidth',1,'color',[0,0,0.6])
    plot(Results(4,:,4),'o','LineWidth',1,'color',[0.5,0,0.6])
    plot(Results(5,:,4),'o','LineWidth',1,'color',[0.6,0,0])
    plot(Results(7,:,4),'o','LineWidth',1,'color',[0.6,0.6,0])
    plot(Results(10,:,4),'o','LineWidth',1,'color',[0,0.6,0])
    plot(Results(14,:,4),'o','LineWidth',1,'color',[0,0.5,0.8])
    xlabel('Active Days','fontsize',12)
    ylabel('Cells (10^{11})','fontsize',12)
    % ylim([0 0.8])
    legend('Continuous','2-day period','3-day period','4-day period'...
        ,'5-day period','7-day period','10-day period','14-day period');
    title('Comparison of normal cells'' minima', 'fontsize', 12)
    
    saveas(fig4, strcat(path,'comp_N_min'), 'fig');
    print(fig4,'-dpng', strcat(path,'comp_N_min.png'));


    fig5=figure(5);
    set(gca,'FontSize',11)
    plot(Results(1,:,5),'ro','LineWidth',1)
    hold on;
    plot(Results(2,:,5),'o','LineWidth',1,'color',[0,0,0]) 
    plot(Results(3,:,5),'o','LineWidth',1,'color',[0,0,0.6])
    plot(Results(4,:,5),'o','LineWidth',1,'color',[0.5,0,0.6])
    plot(Results(5,:,5),'o','LineWidth',1,'color',[0.6,0,0])
    plot(Results(7,:,5),'o','LineWidth',1,'color',[0.6,0.6,0])
    plot(Results(10,:,5),'o','LineWidth',1,'color',[0,0.6,0])
    plot(Results(14,:,5),'o','LineWidth',1,'color',[0,0.5,0.8])
    xlabel('Active Days','fontsize',12)
    ylabel('Cells (10^{11})','fontsize',12)
    % ylim([0 0.35])
    legend('Continuous','2-day period','3-day period','4-day period'...
        ,'5-day period','7-day period','10-day period','14-day period');
    title('Comparison of maximum tumor size', 'fontsize', 12)
    
    saveas(fig5, strcat(path,'comp_T_max'), 'fig');
    print(fig5,'-dpng', strcat(path,'comp_T_max.png'));
   
    
end

