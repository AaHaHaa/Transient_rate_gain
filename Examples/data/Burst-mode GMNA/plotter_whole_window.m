% Plot all pulses
% Note the potential plotting gap between time windows
figure;
yyaxis right;
this_t = 0;
for i = 1:2
    if mod(i,2) == 1
        t_tmp = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'*prop_output(i).dt;
        if i == 1
            this_t = t_tmp - t_tmp(1);
        else
            this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        end
        plot(this_t/1e6,zeros(gain_rate_eqn.woPulse.Nt,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',4);
    else
        t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
        last_t_end = this_t(end);
        this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        plot([last_t_end;this_t]/1e6,abs([0;prop_output(i).fields.forward(:,:,end)]).^2/1e3,'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',4);
    end
    hold on;
end
hold off;
xlabel('Time (µs)')
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);

yyaxis left;
this_t = 0;
for i = 1:2
    t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
    last_t_end = this_t(end);
    if i == 1
        this_t = t_tmp - t_tmp(1);
    else
        this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
    end
    if mod(i,2) == 1
        plot(this_t/1e6,prop_output(i).population(:,:,end)*100,'Color','b','LineStyle','-','Marker','None','linewidth',4);
    else
        plot([last_t_end;this_t]/1e6,[prop_output(i-1).population(end,:,end);prop_output(i).population(:,:,end)]*100,'Color','b','LineStyle','-','Marker','None','linewidth',4);
    end
    hold on;
end
hold off;
%xlabel('Time (\mus)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',30);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
print(gcf,'whole_window.pdf','-dpdf','-bestfit');