close all;
clear i j period active_days Cells_out Drug_out;

save = 1;
max_period = 14;
max_consecutive_days = 7; % from the experiments, no more consecutive days are required

% Initialize only in the first run !!
% 1 - total drug given
% 2 - maximum drug concetration
% 3 - day of tumor eradication
% 4 - minimum normal cells
% 5 - maximum tumor cells
if ~(exist('Results','var'))
    Results=NaN(max_period,max_consecutive_days,5);
end

simTime = 60;
step = 0.1; % double check simulink
sim('model_Banks',simTime);

case_sel=Case_out.data(1); % {1,2,3,4}

% Obtain period/active days
j=0;
for i=1:max_period+1
    if (Pulse_out.data(i)==0)
        j = j+1; % inactive days
        if (Pulse_out.data(i+1)==1) % new period
            break;
        end
    end
end
if (i==15 && j==0) % continuous
    period=1;
    active_days=1;
else
    period = i;
    active_days = i-j;
end

% Drug_out.data = Drug_out.data(Drug_out.data>=0); % remove any possible negative values

% Total drug administered
total_drug = sum(Drug_out.data);
fprintf('[Case %d] [%d/%d] Total drug : %g mg/m^2\n', case_sel, period, active_days, total_drug);

% Maximum drug concentration
max_drug_conc = max(Cells_out.data(:,4));
fprintf('[Case %d] [%d/%d] Max drug concentration : %g mg/L\n', case_sel, period, active_days, max_drug_conc);

% Day of tumor eradication
t_zero = -1;
for i=1:simTime/step
    if Cells_out.data(i,2)<0.005
        t_zero = i*step;
        break;
    end
end
fprintf('[Case %d] [%d/%d] Day of eradication : %g \n', case_sel, period, active_days, t_zero);
% if t_zero == -1 -> ineffective regimen, the tumor wasn't eradicated

Nmin = min(Cells_out.data(:,1));
fprintf('[Case %d] [%d/%d] Minimum normal cells : %g \n', case_sel, period, active_days, Nmin);

Tmax = max(Cells_out.data(:,2));
% fprintf('[Case %d] [%d/%d] Maximum tumor cells : %g \n', case_sel, period, active_days, Tmax);

Results(period,active_days,1) = total_drug;
Results(period,active_days,2) = max_drug_conc;
Results(period,active_days,3) = t_zero;
Results(period,active_days,4) = Nmin;
Results(period,active_days,5) = Tmax;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                               Print                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

path = sprintf('figures\\case %d\\[%d-%d]', case_sel, period, active_days);
%%% Cells figure
fig_cells = figure(1);
% plot(Cells_out.time,Cells_out.data,'LineWidth',1)
if period==1 && active_days==1 
    stairs(Drug_out.time,Drug_out.data,'k','LineWidth',1,'color',[0.35,0.35,0.35])
else
    stairs(Drug_out.time,Drug_out.data,'k','LineWidth',1,'color',[0.4,0.4,0.4])
end
hold on
plot(Cells_out.time,Cells_out.data,'LineWidth',1)


if period==1 && active_days==1 
    title(sprintf('Cell Populations, Drug concentration and Drug input'),'fontsize',12);
else
    title(sprintf('Cell Populations, Drug concentration and Drug input [%d/%d]', period, active_days),'fontsize',12)
end
set(gca,'FontSize',11)
xlabel('Days','fontsize',12)
ylabel('Cells (10^{11}), Drug conc.(mg/L), Drug input(mg/m^2)','fontsize',12)
legend('N','T','I','M','v')

if save == 1 
    saveas(fig_cells, path, 'fig');
    print(fig_cells,'-dpng',strcat(path,'.png'));
end

%%% Drug figure

fig_drug = figure(2);
stairs(Drug_out.time,Drug_out.data,'k','LineWidth',1)

title(sprintf('Drug input. v_{total}: %g', total_drug),'fontsize',12)

set(gca,'FontSize',11)
xlabel('Days','fontsize',12)
ylabel('Drug (mg/m^2)','fontsize',12)
if (case_sel == 4)
    ylim([0 1.8])
end

if save == 1
    saveas(fig_drug, strcat(path,'-d'), 'fig');
    print(fig_drug,'-dpng', strcat(path,'-d.png'));
end

%     plot_periodic_results(case_sel, Results);
