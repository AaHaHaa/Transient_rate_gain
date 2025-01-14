clearvars; close all;

addpath('../../../user_helpers/');

%%
num_pulses = 4;
load(sprintf('GMNA_woASE_%upulses.mat',num_pulses));

figure;
yyaxis left;
confplot(grating_spacings*1e3,I_mean/1e6,I_std/1e6,'Color','b','linewidth',2,'LineStyle','-');
ylabel('Peak power (MW)');
ylim([0,1.46]);
set(gca,'YColor','b');
yyaxis right;
confplot(grating_spacings*1e3,FWHM_mean*1e3,FWHM_std*1e3,'linewidth',2,'LineStyle','-');
ylabel('Pulse duration (fs)');
ylim([20,100]);
xlabel('Grating spacing (mm)');
xlim([min(grating_spacings),max(grating_spacings)]*1e3);
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.08,0.187,0.84,0.73])
print(gcf,sprintf('dechirp_info_%upulses.pdf',num_pulses),'-dpdf','-bestfit');

for ii = [1.04,1.073,1.103]*1e-3
    dechirped_field = pulse_compressor_single( 'Treacy-t',ii,pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    figure;
    plot(t,abs(dechirped_field).^2/1e6,'linewidth',4,'Color','b');
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-15,15]);
    ylim([0,1.3]);
    print(gcf,sprintf('DAout_%4.2fmm_%upulses.pdf',ii*1e3,num_pulses),'-dpdf');
end

%%
%{
num_pulses = 100;
load(sprintf('GMNA_woASE_%upulses.mat',num_pulses));

figure;
yyaxis left;
confplot(grating_spacings*1e3,I_mean/1e6,I_std/1e6,'Color','b','linewidth',2,'LineStyle','-');
ylabel('Peak power (MW)');
ylim([0,1.46]);
set(gca,'YColor','b');
yyaxis right;
confplot(grating_spacings*1e3,FWHM_mean*1e3,FWHM_std*1e3,'linewidth',2,'LineStyle','-');
ylabel('Pulse duration (fs)');
ylim([20,100]);
xlabel('Grating spacing (mm)');
xlim([min(grating_spacings),max(grating_spacings)]*1e3);
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.08,0.187,0.84,0.73])
print(gcf,sprintf('dechirp_info_%upulses.pdf',num_pulses),'-dpdf','-bestfit');

for ii = [1.65,1.71,1.7275,1.75,1.8]*1e-3
    dechirped_field = pulse_compressor_single( 'Treacy-t',ii,pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    figure;
    plot(t,abs(dechirped_field).^2/1e6,'linewidth',4,'Color','b');
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-375,375]);
    ylim([0,1.3]);
    print(gcf,sprintf('DAout_%4.2fmm_%upulses.pdf',ii*1e3,num_pulses),'-dpdf');
end
%}
%%
num_pulses = 200;
load(sprintf('GMNA_woASE_%upulses.mat',num_pulses));

figure;
yyaxis left;
confplot(grating_spacings*1e3,I_mean/1e6,I_std/1e6,'Color','b','linewidth',2,'LineStyle','-');
ylabel('Peak power (MW)');
ylim([0,1.46]);
set(gca,'YColor','b');
yyaxis right;
confplot(grating_spacings*1e3,FWHM_mean*1e3,FWHM_std*1e3,'linewidth',2,'LineStyle','-');
ylabel('Pulse duration (fs)');
ylim([20,100]);
xlabel('Grating spacing (mm)');
xlim([min(grating_spacings),max(grating_spacings)]*1e3);
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.08,0.187,0.84,0.73])
print(gcf,sprintf('dechirp_info_%upulses.pdf',num_pulses),'-dpdf','-bestfit');

for ii = [1.043,1.058,1.079,1.124,1.151,1.196]*1e-3
    dechirped_field = pulse_compressor_single( 'Treacy-t',ii,pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    figure;
    plot(t,abs(dechirped_field).^2/1e6,'linewidth',4,'Color','b');
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-800,800]);
    ylim([0,1.6]);
    print(gcf,sprintf('DAout_%4.2fmm_%upulses.pdf',ii*1e3,num_pulses),'-dpdf');
end