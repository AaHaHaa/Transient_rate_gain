close all; clearvars;

%% Simulation data
load('singleCPA_Yilmaz_1000kHz.mat');
sim_Nt_10pulses = Nt;
sim_dt_10pulses = [prop_output.dt];
sim_rep_rate_10pulses = burst_rep_rate;
sim_population_10pulses = [prop_output.population]; % this symbol puts data of different windows into the 2nd dimension
sim_t_10pulses = t/1e3; % ns
sim_I_10pulses = abs(prop_output(2).fields.forward(:,:,end)).^2/1e3; sim_I_10pulses = sim_I_10pulses/max(sim_I_10pulses);

load('singleCPA_Yilmaz_500kHz.mat');
sim_t_20pulses = t/1e3; % ns
sim_I_20pulses = abs(prop_output(2).fields.forward(:,:,end)).^2/1e3; sim_I_20pulses = sim_I_20pulses/max(sim_I_20pulses);

load('singleCPA_Yilmaz_200kHz.mat');
sim_t_50pulses = t/1e3; % ns
sim_I_50pulses = abs(prop_output(2).fields.forward(:,:,end)).^2/1e3; sim_I_50pulses = sim_I_50pulses/max(sim_I_50pulses);

%% Paper data
paper_10pulses = readmatrix('1 MHz 10 pulses.csv');
paper_t_10pulses = paper_10pulses(:,1) - 80; % ns
paper_I_10pulses = paper_10pulses(:,2); % norm.

[paper_t_10pulses,ia] = unique(paper_t_10pulses);
paper_I_10pulses = paper_I_10pulses(ia);

background = mean(paper_I_10pulses(end-7000:end));
paper_I_10pulses = paper_I_10pulses-background; paper_I_10pulses = paper_I_10pulses/max(paper_I_10pulses);

paper_20pulses = readmatrix('500 kHz 20 pulses.csv');
paper_t_20pulses = paper_20pulses(:,1) - 250; % ns
paper_I_20pulses = paper_20pulses(:,2); % norm.

[paper_t_20pulses,ia] = unique(paper_t_20pulses);
paper_I_20pulses = paper_I_20pulses(ia);

background = mean(paper_I_20pulses(end-3000:end));
paper_I_20pulses = paper_I_20pulses-background; paper_I_20pulses = paper_I_20pulses/max(paper_I_20pulses);

paper_50pulses = readmatrix('200 kHz 50 pulses.csv');
paper_t_50pulses = paper_50pulses(:,1) - 400; % ns
paper_I_50pulses = paper_50pulses(:,2); % norm.

[paper_t_50pulses,ia] = unique(paper_t_50pulses);
paper_I_50pulses = paper_I_50pulses(ia);

background = mean(paper_I_50pulses(end-1500:end));
paper_I_50pulses = paper_I_50pulses-background; paper_I_50pulses = paper_I_50pulses/max(paper_I_50pulses);

%% Shift the paper data temporally
sim_e10 = envelope(sim_I_10pulses,10,'peak');
paper_e10 = envelope(paper_I_10pulses,10,'peak');
opt_fun_10pulses = @(dt) overlap(dt,sim_t_10pulses,sim_e10,paper_t_10pulses,paper_e10);
opt_t_10pulses = fminsearch(opt_fun_10pulses,0);
paper_t_10pulses = paper_t_10pulses + opt_t_10pulses;

sim_e20 = envelope(sim_I_20pulses,10,'peak');
paper_e20 = envelope(paper_I_20pulses,10,'peak');
opt_fun_20pulses = @(dt) overlap(dt,sim_t_20pulses,sim_e20,paper_t_20pulses,paper_e20);
opt_t_20pulses = fminsearch(opt_fun_20pulses,0);
paper_t_20pulses = paper_t_20pulses + opt_t_20pulses;

sim_e50 = envelope(sim_I_50pulses,10,'peak');
paper_e50 = envelope(paper_I_50pulses,10,'peak');
opt_fun_50pulses = @(dt) overlap(dt,sim_t_50pulses,sim_e50,paper_t_50pulses,paper_e50);
opt_t_50pulses = fminsearch(opt_fun_50pulses,0);
paper_t_50pulses = paper_t_50pulses + opt_t_50pulses;

%% Plot each
figure;
plot(paper_t_10pulses,paper_I_10pulses,'linewidth',3,'Color','r');
hold on;
plot(sim_t_10pulses,sim_I_10pulses,'linewidth',6,'Color','b');
hold off;
xlabel('Time (ns)');
ylabel('Power (norm.)');
xlim([-75,75]);
ylim([-0.05,1]);
set(gca,'fontsize',30);
print(gcf,'1 MHz.pdf','-dpdf');

figure;
plot(paper_t_20pulses,paper_I_20pulses,'linewidth',2,'Color','r');
hold on;
plot(sim_t_20pulses,sim_I_20pulses,'linewidth',4,'Color','b');
hold off;
xlabel('Time (ns)');
ylabel('Power (norm.)');
xlim([-150,150]);
ylim([-0.05,1]);
set(gca,'fontsize',30);
print(gcf,'500 kHz.pdf','-dpdf');

figure;
plot(paper_t_50pulses,paper_I_50pulses,'linewidth',1.5,'Color','r');
hold on;
plot(sim_t_50pulses,sim_I_50pulses,'linewidth',2,'Color','b');
hold off;
xlabel('Time (ns)');
ylabel('Power (norm.)');
xlim([-375,375]);
ylim([-0.05,1]);
set(gca,'fontsize',30);
print(gcf,'200 kHz.pdf','-dpdf');

%% Plot all pulses for the 1-MHz case
figure;
yyaxis right;
this_t = 0;
for i = 1:2
    if mod(i,2) == 1
        t_tmp = (-sim_Nt_10pulses/2:sim_Nt_10pulses/2-1)'*sim_dt_10pulses(i);
        if i == 1
            this_t = t_tmp - t_tmp(1);
        else
            this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        end
        plot(this_t/1e6,zeros(sim_Nt_10pulses,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    else
        t_tmp = (-size(sim_population_10pulses,1)/2:size(sim_population_10pulses,1)/2-1)'*sim_dt_10pulses(i);
        last_t_end = this_t(end);
        this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        plot([last_t_end;this_t]/1e6,[0;sim_I_10pulses],'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
xlabel('Time (µs)')
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);

yyaxis left;
this_t = 0;
for i = 1:2
    t_tmp = (-size(sim_population_10pulses,1)/2:size(sim_population_10pulses,1)/2-1)'*sim_dt_10pulses(i);
    last_t_end = this_t(end);
    if i == 1
        this_t = t_tmp - t_tmp(1);
    else
        this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
    end
    if mod(i,2) == 1
        plot(this_t/1e6,sim_population_10pulses(:,i,end)*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
    else
        plot([last_t_end;this_t]/1e6,[sim_population_10pulses(end,i-1,end);sim_population_10pulses(:,i,end)]*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
%xlabel('Time (\mus)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
xlim([0,1/sim_rep_rate_10pulses*1e6]);
print(gcf,'all_pulses.pdf','-dpdf','-bestfit','-vector');

function result = overlap(t,t1,I1,t2,I2)

I2 = interp1(t2+t,I2,t1,'linear',0);

result = -trapz(t1,I1.*I2);

end