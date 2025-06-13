clearvars; close all;

num_pulses = 200;

load(sprintf('GMNA_woASE_%upulses_preshaping30.mat',num_pulses));

figure;
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',3,'Color','b');
ylim([4.7,5.2]);
ylabel('N_1 (%)');
set(gca,'YColor','b');
yyaxis right;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',0.5);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
xlim([-5,5]*num_pulses);

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.085,0.187,0.845,0.73])

print(gcf,sprintf('Aout_%upulses_preshaping.pdf',num_pulses),'-dpdf','-bestfit');

figure;
plot(t,abs(prop_output(2).fields.forward(:,:,1)).^2,'linewidth',0.5,'Color','k');
xlabel('Time (ps)');
set(gca,'fontsize',40,'YTick',[]);
xlim([-5,5]*num_pulses);
print(gcf,sprintf('Ain_%upulses_preshaping.pdf',num_pulses),'-dpdf');