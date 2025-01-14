close all; clearvars;

grating_angles = 19;%19:21;

%% non-preshaping
load(sprintf('GMNA_woASE_%upulses.mat',200));

num_save = length(prop_output(1).z);

each_energy{1} = zeros(num_save,num_pulses);
each_peak_power{1} = zeros(num_save,num_pulses);
each_nonlinear_phase{1} = zeros(1,num_pulses);
for i = 1:num_save
    [~,pulse_span] = findpeaks(-abs(prop_output(2).fields.forward(:,:,i)).^2,'MinPeakDistance',(DT-3)/dt,'MinPeakHeight',-100,'MinPeakWidth',1/dt,'MinPeakProminence',100);
    pulse_span = [1;pulse_span;Nt];
    for j = 1:num_pulses
        idx = pulse_span(j):pulse_span(j+1);
        each_energy{1}(i,j) = trapz(t(idx),abs(prop_output(2).fields.forward(idx,:,i)).^2)/1e3;
        each_peak_power{1}(i,j) = max(abs(prop_output(2).fields.forward(idx,:,i)).^2);
    end
end
% Dechirp based on the preset grating angle
for jj = 1:length(grating_angles)
    for ii = 1:length(grating_spacings)
        dechirped_field_i = pulse_compressor_single( 'Treacy-t',grating_spacings(ii),pi/180*grating_angles(jj),sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    
        I = abs(dechirped_field_i).^2;
        threshold = max(I)/3;
        [this_max_I,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2,'MinPeakDistance',DT*0.9);
    
        I_mean(ii) = mean(this_max_I);
        I_std(ii) = std(this_max_I);
    end
    variation_idx = I_std./I_mean;
    variation{jj,1} = min(variation_idx(I_mean > max(I_mean)*0.5));
end
c = 299792458;
nonlin_const = fiber.n2*(sim.f0*2*pi*1e12)/c*fiber.SR; % W^-1 m
each_nonlinear_phase{1} = nonlin_const*trapz(prop_output(1).z,each_peak_power{1},1);

%% preshaping
ratio = [1,1.5,1.6,1.7,1.8,1.9,2,2.5,2.8,3,3.3,3.5,4];

for k = 1:length(ratio)-1
    load(sprintf('GMNA_woASE_%upulses_preshaping%u.mat',200,ratio(k+1)*10));

    each_energy{k+1} = zeros(num_save,num_pulses);
    each_peak_power{k+1} = zeros(num_save,num_pulses);
    each_nonlinear_phase{k+1} = zeros(1,num_pulses);
    for i = 1:num_save
        [~,pulse_span] = findpeaks(-abs(prop_output(2).fields.forward(:,:,i)).^2,'MinPeakDistance',(DT-3)/dt,'MinPeakHeight',-100,'MinPeakWidth',1/dt,'MinPeakProminence',100);
        pulse_span = [1;pulse_span;Nt];
        for j = 1:num_pulses
            idx = pulse_span(j):pulse_span(j+1);
            each_energy{k+1}(i,j) = trapz(t(idx),abs(prop_output(2).fields.forward(idx,:,i)).^2)/1e3;
            each_peak_power{k+1}(i,j) = max(abs(prop_output(2).fields.forward(idx,:,i)).^2);
        end
    end
    each_nonlinear_phase{k+1} = nonlin_const*trapz(prop_output(1).z,each_peak_power{k+1},1);

    % Dechirp based on the preset grating angle
    for jj = 1:length(grating_angles)
        for ii = 1:length(grating_spacings)
            dechirped_field_i = pulse_compressor_single( 'Treacy-t',grating_spacings(ii),pi/180*grating_angles(jj),sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
        
            I = abs(dechirped_field_i).^2;
            threshold = max(I)/3;
            [this_max_I,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2,'MinPeakDistance',DT*0.9);
        
            I_mean(ii) = mean(this_max_I);
            I_std(ii) = std(this_max_I);
        end
        variation_idx = I_std./I_mean;
        variation{jj,k+1} = min(variation_idx(I_mean > max(I_mean)*0.5));
    end
end

variation = cell2mat(variation);

%% Overall analysis
ccc = distinguishable_colors(4);
figure;
plot(each_energy{ratio==1  }(end,:),'linewidth',2,'Color',ccc(1,:));
hold on;
plot(each_energy{ratio==1.7}(end,:),'linewidth',2,'Color',ccc(2,:));
plot(each_energy{ratio==2  }(end,:),'linewidth',2,'Color',ccc(3,:));
plot(each_energy{ratio==3  }(end,:),'linewidth',2,'Color',ccc(4,:));
ylabel('Energy (nJ)');
xlabel('Pulse');
set(gca,'fontsize',20);
print(gcf,'energy_pulse.pdf','-dpdf');

figure;
plot(each_nonlinear_phase{ratio==1  }(end,:),'linewidth',2,'Color',ccc(1,:));
hold on;
plot(each_nonlinear_phase{ratio==1.7}(end,:),'linewidth',2,'Color',ccc(2,:));
plot(each_nonlinear_phase{ratio==2  }(end,:),'linewidth',2,'Color',ccc(3,:));
plot(each_nonlinear_phase{ratio==3  }(end,:),'linewidth',2,'Color',ccc(4,:));
ylabel('Nonlinear phase (rad)');
xlabel('Pulse');
set(gca,'fontsize',20);
print(gcf,'nonlinear_phase_pulse.pdf','-dpdf');

figure;
ccc = distinguishable_colors(length(grating_angles));
h = plot(ratio,variation.','linewidth',2);
set(h(1),'Color',ccc(1,:));
xlabel('Seed''s preshaped energy ratio');
ylabel('Peak-power variation');
set(gca,'fontsize',20);
YLim = get(gca,'YLim'); ylim([min(variation(:)),YLim(2)]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.085,0.187,0.845,0.73]);
if length(h) > 1
    for jj = 2:size(variation,1)
        set(h(jj),'Color',ccc(jj,:));
    end
    legend(arrayfun(@(x) [num2str(x),char(176)],grating_angles,'UniformOutput',false));
    print(gcf,'preshaped_study.pdf','-dpdf','-bestfit');
else
    print(gcf,sprintf('preshaped_study_GratingAngle%u.pdf',grating_angles),'-dpdf','-bestfit');
end

%% Detailed (each) analysis
ratio_to_plot = [1,1.7,3,4];
for k = 1:length(ratio_to_plot)
    if ratio_to_plot(k) == 1
        load(sprintf('GMNA_woASE_%upulses.mat',200));
    else
        load(sprintf('GMNA_woASE_%upulses_preshaping%u.mat',200,ratio_to_plot(k)*10));
    end

    figure;
    plot(t,abs(prop_output(2).fields.forward(:,:,1)).^2/1e3,'linewidth',0.5,'Color','k');
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-5,5]*num_pulses);
    ylim([0,0.76]);
    print(gcf,sprintf('Ain_%upulses_preshaping%u.pdf',num_pulses,ratio_to_plot(k)*10),'-dpdf');

    figure;
    plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',0.5,'Color',[0.8510,0.3255,0.0980]);
    xlabel('Time (ps)');
    set(gca,'fontsize',40,'YTick',[]);
    xlim([-5,5]*num_pulses);
    ylim([0,16.5]);
    print(gcf,sprintf('Aout_%upulses_preshaping%u.pdf',num_pulses,ratio_to_plot(k)*10),'-dpdf');
end