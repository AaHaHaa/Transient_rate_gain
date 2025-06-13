
clearvars; close all;

addpath('../../../user_helpers/');

DT_all = 0.1:0.1:7;

initial_guess = 0;
max_I = zeros(length(DT_all),1);
energy_all = zeros(length(DT_all),1);
for iii = length(DT_all):-1:1
    DT = DT_all(iii);
    
    load(sprintf('GMNA_woASE_%2.1fps_2pulses.mat',DT));
    
    [grating_spacing,~,dechirped_field] = pulse_compressor_bursts( 'Treacy-t',pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6,initial_guess,false,true );
    max_I(iii) = max(abs(dechirped_field).^2);
    
    initial_guess = grating_spacing;
    energy_all(iii) = energy(end);
end

figure;
yyaxis left;
plot(DT_all',max_I/1e6,'linewidth',2,'Color','b');
set(gca,'fontsize',20,'YColor','b');
xlabel('Seed-pulse separation (ps)');
ylabel('Peak power (MW)');
yyaxis right;
plot(DT_all',energy_all,'linewidth',2);
YLim = get(gca,'YLim');
ylim([70,YLim(2)]);
ylabel('Energy (nJ)');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.08,0.187,0.84,0.73])
print(gcf,'peak_power.pdf','-dpdf','-bestfit');