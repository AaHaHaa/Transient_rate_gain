clearvars; close all;

addpath('../../../user_helpers/');

grating_angle = 21;
preshaped_ratio_all = [1.7,2,3];

num_pulses = 200;

figure;
ccc = distinguishable_colors(length(preshaped_ratio_all));
for ii = 1:length(preshaped_ratio_all)
    load(sprintf('GMNA_woASE_%upulses_preshaping%u.mat',num_pulses,preshaped_ratio_all(ii)*10));

    [grating_spacing,~,dechirped_field] = pulse_compressor_bursts( 'Treacy-t',pi/180*grating_angle,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );

    % Dechirp based on the preset grating angle
    for i = 1:length(grating_spacings)
        dechirped_field_i = pulse_compressor_single( 'Treacy-t',grating_spacings(i),pi/180*grating_angle,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
    
        I = abs(dechirped_field_i).^2;
        threshold = max(I)/3;
        [this_max_I,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2,'MinPeakDistance',DT*0.9);
    
        I_mean(i) = mean(this_max_I);
        I_std(i) = std(this_max_I);
    end

    yyaxis left;
    hhh1(ii) = plot(grating_spacings*1e3,I_std,'linewidth',2,'Color',ccc(ii,:),'linestyle','-','MarkerEdgeColor','None'); hold on;
    
    idx = find(grating_spacings>grating_spacing,1);
    plot(grating_spacings(idx)*1e3,I_std(idx),'MarkerFaceColor',ccc(ii,:),'MarkerSize',20,'Marker','o','MarkerEdgeColor','k','LineStyle','None');

    yyaxis right;
    hhh2(ii) = plot(grating_spacings*1e3,I_mean/1e6,'linewidth',2,'Color',ccc(ii,:),'linestyle','--'); hold on;
end
min_g = min(grating_spacings*1e3);
max_g = max(grating_spacings*1e3);
yyaxis left;
hold off;
set(gca,'fontsize',20);
xlabel('Grating spacing (mm)');
ylabel('Peak-power variation');
set(gca,'YColor','b');
%ylabel('\sigma_{peak power} (MW)');
xlim([min_g,max_g]);
%ylim([0,0.4]);
ylim([0,4e5]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
set(gca,'Position',[0.085,0.187,0.83,0.73]);
l = legend(hhh1,'1.7','2','3'); set(l,'location','northwest');
yyaxis right;
ylabel('Peak power (MW)');
print(gcf,sprintf('std(I)_preshaping%u_GratingAngle%u.pdf',preshaped_ratio_all(ii)*10,grating_angle),'-dpdf','-bestfit');