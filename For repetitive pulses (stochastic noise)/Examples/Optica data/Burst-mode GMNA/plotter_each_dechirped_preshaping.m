clearvars; close all;

%%
num_pulses = 200;
load(sprintf('GMNA_woASE_%upulses_preshaping30.mat',num_pulses));

figure;
yyaxis left;
confplot(grating_spacings*1e3,I_mean/1e6,I_std/1e6,'Color','b','linewidth',2,'LineStyle','-');
ylabel('Peak power (MW)');
ylim([0,1.46]);
set(gca,'YColor','b');
yyaxis right;
confplot(grating_spacings*1e3,FWHM_mean*1e3,FWHM_std*1e3,'linewidth',2,'LineStyle','-');
ylabel('Pulse duration (fs)');
ylim([20,180]);
xlabel('Grating spacing (mm)');
xlim([min(grating_spacings),max(grating_spacings)]*1e3);
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.08,0.187,0.84,0.73])
print(gcf,sprintf('dechirp_info_%upulses_preshaping.pdf',num_pulses),'-dpdf','-bestfit');

for ii = [1.106,1.121,1.139,1.163,1.175]*1e-3
    dechirped_field = pulse_compressor_single( 'Treacy-t',ii,pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    figure;
    plot(t,abs(dechirped_field).^2/1e6,'linewidth',4,'Color','b');
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-800,800]);
    ylim([0,1.6]);
    print(gcf,sprintf('DAout_%4.2fmm_%upulses_preshaping.pdf',ii*1e3,num_pulses),'-dpdf');
end