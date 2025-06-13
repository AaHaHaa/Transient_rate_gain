clearvars; close all;

addpath('../../../user_helpers/');

each(7);
each(3);
each(1);
each(0.5);

%%
function each(DT)
    load(sprintf('GMNA_woASE_%2.1fps_2pulses.mat',DT));

    figure;
    plot(t,abs(prop_output(2).fields.forward(:,:,1)).^2,'linewidth',10,'Color','k');
    set(gca,'fontsize',40,'YTick',[],'XTick',[],'Color','None','XColor','None','YColor','None');
    xlim([-4.5,4.5]);
    ylim([0,600]);
    box off;
    print(gcf,sprintf('Ain_%2.1fps_2pulses.pdf',DT),'-dpdf');
end