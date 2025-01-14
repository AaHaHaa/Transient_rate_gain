close all; clearvars;

%% non-preshaping
load(sprintf('GMNA_woASE_%upulses.mat',200));

factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

each_spectrum{1} = zeros(Nt,num_pulses);
[~,pulse_span] = findpeaks(-abs(prop_output(2).fields.forward(:,:,end)).^2,'MinPeakDistance',(DT-3)/dt,'MinPeakHeight',-100,'MinPeakWidth',1/dt,'MinPeakProminence',100);
pulse_span = [1;pulse_span;Nt];
for j = 1:num_pulses
    idx = pulse_span(j):pulse_span(j+1);
    pulse_field = zeros(Nt,1);
    pulse_field(idx) = prop_output(2).fields.forward(idx,:,end);
    each_spectrum{1}(:,j) = abs(fftshift(ifft(pulse_field),1)).^2*factor_correct_unit.*factor;
end

%% preshaping
ratio = [1,2,3,4];

for k = 1:length(ratio)-1
    load(sprintf('GMNA_woASE_%upulses_preshaping%u.mat',200,ratio(k+1)*10));

    each_spectrum{k+1} = zeros(Nt,num_pulses);
    [~,pulse_span] = findpeaks(-abs(prop_output(2).fields.forward(:,:,end)).^2,'MinPeakDistance',(DT-3)/dt,'MinPeakHeight',-100,'MinPeakWidth',1/dt,'MinPeakProminence',100);
    pulse_span = [1;pulse_span;Nt];
    for j = 1:num_pulses
        idx = pulse_span(j):pulse_span(j+1);
        pulse_field = zeros(Nt,1);
        pulse_field(idx) = prop_output(2).fields.forward(idx,:,end);
        each_spectrum{k+1}(:,j) = abs(fftshift(ifft(pulse_field),1)).^2*factor_correct_unit.*factor;
    end
end

%% Plot
for k = 1:length(ratio)
    figure;
    plot(lambda,each_spectrum{k},'linewidth',2,'Color',[0.5,0.5,0.5]);
    hold on;
    plot(lambda,each_spectrum{k}(:,1),'linewidth',2,'Color','b');
    plot(lambda,each_spectrum{k}(:,end),'linewidth',2,'Color','r');
    hold off;
    set(gca,'fontsize',20);
    xlim([950,1200]);
    xlabel('Wavelength (nm)');
    ylabel('PSD (nJ/nm)');
    print(gcf,sprintf('spectral_variations_preshaping%u.pdf',ratio(k)*10),'-dpdf');
end