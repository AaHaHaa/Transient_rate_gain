close all; clearvars;

total_rt = 3000;
rt = (1:total_rt)';

energy_all = zeros(total_rt,1);
peak_power_all = zeros(total_rt,1);
population1_all = zeros(total_rt,1);
population2_all = zeros(total_rt,1);
population3_all = zeros(total_rt,1);
population4_all = zeros(total_rt,1);

name_idx = 1000:1000:total_rt;
for i = 1:length(name_idx)
    load(sprintf('oscillator_%u.mat',name_idx(i)));

    energy_all(1000*(i-1)+1:1000*i) = sum(output_energy,2);
    peak_power_all(1000*(i-1)+1:1000*i) = max([squeeze(max(abs(output_field{1}).^2,[],1)),squeeze(max(abs(output_field{2}).^2,[],1))],[],2);
    population1_all(1000*(i-1)+1:1000*i) = population1;
    population2_all(1000*(i-1)+1:1000*i) = population2;
    population3_all(1000*(i-1)+1:1000*i) = population3;
    population4_all(1000*(i-1)+1:1000*i) = population4;
end

fig = figure;
ax1 = axes();
ax2 = copyobj(ax1,fig); % create two identical overlaying axes
ax3 = copyobj(ax1,fig); % create two identical overlaying axes
linkaxes([ax1,ax2,ax3],'x'); % link the x axes
% Energy and peak power
plot(ax1,rt/1e3,energy_all/1e3,'Color','b','linewidth',4); ylabel(ax1,'Energy (\muJ)'); xlim(ax1,[0,5]); ylim(ax1,[0,13]);
plot(ax2,rt/1e3,peak_power_all/1e3,'Color','r','linewidth',4); ylabel(ax2,'Peak power (kW)'); xlim(ax2,[0,5]); ylim(ax2,[0,90]);
% Populations
h = plot(ax3,rt/1e3,[population1_all,population2_all,population3_all,population4_all]*100,'linewidth',4); ylabel(ax3,'N_1 (%)'); xlim(ax3,[0,5]); ylim(ax3,[0,55]); set(ax3, 'YAxisLocation', 'right'); % set y-axes to right side
set(h(1),'Color',0*ones(1,3)); set(h(2),'Color',0.4*ones(1,3)); set(h(3),'Color',0.6*ones(1,3)); set(h(4),'Color',0.7*ones(1,3));
xlabel(ax1,'Roundtrip (x10^3)'); xlabel(ax2,''); xlabel(ax3,'');
set([ax1,ax2,ax3],'Color', 'None','fontsize',20);
set(ax1,'YColor','b'); set(ax2,'YColor','r'); set(ax3,'YColor','k');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3),pos(4)]); % make the figure "wider"
pos = ax1.Position; set([ax1,ax2,ax3], 'Position', [0.143,pos(2),0.77,pos(4)]); % make the axes 10% more narrow
% Extend the ax1 y-axes ticklabels rightward a bit by adding space
ax1.YTickLabel = strcat(ax1.YTickLabel,{'            '});
exportgraphics(gcf, 'evolution.pdf', 'ContentType', 'vector'); % print command doesn't work for figures with a far-away ylabel