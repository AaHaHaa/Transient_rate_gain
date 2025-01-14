clearvars; close all;

addpath('../../../user_helpers/');

x0 = each(7,[-4.5,4.5],0);
x0 = each(6,[-4,4],x0);
x0 = each(5,[-3.5,3.5],x0);
x0 = each(4,[-3,3],x0);
x0 = each(3,[-2.5,2.5],x0);
x0 = each(2,[-2,2],x0);
x0 = each(1,[-0.2,0.2],x0);
x0 = each(0.1,[-0.2,0.2],x0);

%%
function grating_spacing = each(DT,Xlim,initial_guess)
    load(sprintf('GMNA_woASE_%2.1fps_2pulses.mat',DT));

    figure;
    yyaxis left;
    plot(t,prop_output(2).population(:,:,end)*100,'linewidth',3,'Color','b');
    ylabel('N_1 (%)');
    set(gca,'YColor','b');
    yyaxis right;
    plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',3);
    xlabel('Time (ps)'); ylabel('Power (kW)');
    ylim([0,40]);
    set(gca,'fontsize',25);
    xlim([-10,10]);
    print(gcf,sprintf('Aout_%2.1fps_2pulses.pdf',DT),'-dpdf');
    
    set(gca,'Units','pixel');
    pos0 = get(gca,'Position');
    
    figure;
    yyaxis left;
    plot(t,prop_output(2).Power.pump.forward(:,:,end),'linewidth',3,'Color','b');
    ylabel('Pump power (W)');
    set(gca,'YColor','b');
    yyaxis right;
    plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',3);
    xlabel('Time (ps)'); ylabel('Power (kW)');
    ylim([0,40]);
    set(gca,'fontsize',25);
    xlim([-10,10]);
    yyaxis left;
    current_YTick = get(gca,'YTick');
    new_YTick = unique(arrayfun(@(x)str2double(sprintf('%.5f',x)),current_YTick));
    set(gca,'Units','pixel');
    set(gca,'YTick',new_YTick,'Position',pos0);
    print(gcf,sprintf('Pump_%2.1fps_2pulses.pdf',DT),'-dpdf');
    
    spectrum = abs(fftshift(ifft(prop_output(2).fields.forward(:,:,end)),1)).^2;
    factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e3" is to make pJ into nJ
    c = 299792.458;
    factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain
    spectrum = spectrum*factor_correct_unit.*factor;
    
    figure;
    plot(lambda,spectrum,'linewidth',3,'Color','b');
    xlabel('Wavelength (nm)');
    ylabel('PSD (nJ/nm)');
    set(gca,'fontsize',25);
    xlim([950,1350]);
    print(gcf,sprintf('Afout_%2.1fps_2pulses.pdf',DT),'-dpdf');

    [grating_spacing,~,dechirped_field] = pulse_compressor_bursts( 'Treacy-t',pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6,initial_guess );
    [~,max_idx] = max(abs(dechirped_field).^2);
    if DT < 2
        dechirped_field = circshift(dechirped_field,-(max_idx-Nt/2));
    end
    
    figure;
    plot(t,abs(dechirped_field).^2/1e6,'linewidth',3,'Color','k');
    xlabel('Time (ps)');
    ylabel('Power (MW)');
    set(gca,'fontsize',25);
    xlim(Xlim);
    ylim([0,2.5]);
    print(gcf,sprintf('DAout_%2.1fps_2pulses.pdf',DT),'-dpdf');
end