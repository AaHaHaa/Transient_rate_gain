close all; clearvars;

%% Paper data
paper_data = readmatrix('Paschotta_data.csv');
paper_t = paper_data(:,1); % ns
paper_I = paper_data(:,2); % kW

%% Plot with simulation data
load('transient.mat');

t_woPulse = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps
t_woPulse = t_woPulse - t_woPulse(1);
figure;
orange = [0.85,0.325,0.098];
plot(t_woPulse/1e9,prop_output(1).population(:,:,end)*100,'linewidth',2,'Color','b');
xlabel('Time (ms)'); ylabel('N_1 (%)');
set(gca,'fontsize',14,'YColor','b');
pos = get(gcf,'Position');
ylim([4,22]);
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
print(gcf,'transient_window1.pdf','-dpdf','-bestfit');
%exportgraphics(gcf,"transient_window1.pdf",'ContentType','vector')

figure;
yyaxis left;
plot(t/1e3,prop_output(2).population(:,:,end)*100,'linewidth',2,'Color','b');
ylabel('N_1 (%)');
ylim([4,22]);
set(gca,'YColor','b');
yyaxis right;
h1 = plot(t/1e3,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2,'Color',orange);
xlabel('Time (ns)');
ylabel('Power (kW)');

hold on
load('./With steady-state rate-gain code/Steady-state.mat');
h2 = plot(t/1e3,abs(prop_output.fields(:,:,end)).^2/1e3,'linewidth',2,'Color',[0 0.4470 0.7410],'LineStyle','--');

h3 = plot(paper_t,paper_I,'linewidth',2,'Color','r');

xlim([-2,2]);
ylim([0,160]);

l = legend([h1,h2,h3],'Ours','Steady-state','Prior work');
set(l,'fontsize',18);
set(gca,'fontsize',20);
set(gca,'YColor',orange);

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*1.2,pos(4)]);

print(gcf,'transient.pdf','-dpdf');