clearvars; close all;

num_pulses_all = [4,10,50,100,200];

min_g = 1000; % random large number
max_g = 0;
ccc = distinguishable_colors(length(num_pulses_all));
figure;
for ii = 1:length(num_pulses_all)
    load(sprintf('GMNA_woASE_%upulses.mat',num_pulses_all(ii)));
    
    hhh(ii) = plot(grating_spacings*1e3,I_std/1e6,'linewidth',2,'Color',ccc(ii,:)); hold on;
    
    idx = find(grating_spacings>grating_spacing,1);
    plot(grating_spacings(idx)*1e3,I_std(idx)/1e6,'MarkerFaceColor',ccc(ii,:),'MarkerSize',20,'Marker','o','MarkerEdgeColor','k','LineStyle','None');
    
    min_g = min(min(grating_spacings*1e3),min_g);
    max_g = max(max(grating_spacings*1e3),max_g);
end
hold off;
set(gca,'fontsize',20);
xlabel('Grating spacing (mm)');
ylabel('\sigma_{peak power} (MW)');
xlim([min_g,max_g]);
l = legend(hhh,'4 pulses (b)','10 pulses','50 pulses','100 pulses','200 pulses (c)');
set(l,'location','northwest');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.085,0.187,0.88,0.73]);
print(gcf,'std(I).pdf','-dpdf','-bestfit');