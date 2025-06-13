clearvars; close all;

each(10,3);
each(100,1);
each(200,0.5);

%%
function each(num_pulses,linewidth)
    load(sprintf('GMNA_woASE_%upulses.mat',num_pulses));

    figure;
    yyaxis left;
    plot(t,prop_output(2).population(:,:,end)*100,'linewidth',3,'Color','b');
    YLim = get(gca,'YLim'); ylim([YLim(1),YLim(2)*1.001]);
    ylabel('N_1 (%)');
    set(gca,'YColor','b');
    yyaxis right;
    plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',linewidth);
    xlabel('Time (ps)'); ylabel('Power (kW)');
    set(gca,'fontsize',25);
    xlim([-5,5]*num_pulses);
    print(gcf,sprintf('Aout_%upulses.pdf',num_pulses),'-dpdf');

    disp((max(prop_output(2).population(:,:,end))-min(prop_output(2).population(:,:,end)))*100);
end