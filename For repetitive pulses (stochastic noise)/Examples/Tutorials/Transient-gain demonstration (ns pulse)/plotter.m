close all; clearvars;

%% Plot with simulation data
load('transient_gain.mat');

% Plot all pulses
% Note the potential plotting gap between time windows, so plotting below
% requires care to fill it.
fig = figure;
save_point = size(prop_output(2).fields.forward,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for j = 1:save_point                                                                                                                                                                                                                                                                                                                                                    
    figure(fig);
    yyaxis right;
    this_t = 0;
    for i = 1:2
        if mod(i,2) == 1
            t_tmp = (-Nt/2:Nt/2-1)'*prop_output(i).dt;
            if i == 1
                this_t = t_tmp - t_tmp(1);
            else
                this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
            end
            plot(this_t/1e6,zeros(Nt,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
        else
            t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
            last_t_end = this_t(end);
            this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
            plot([last_t_end;this_t]/1e6,[0;abs(prop_output(i).fields.forward(:,:,j)).^2/1e3],'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
        end
        hold on;
    end
    hold off;
    xlabel('Time (µs)')
    ylabel('Power (kW)');
    ylim([0,1]);
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
            plot(this_t/1e6,prop_output(i).population(:,:,j)*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
        else
            plot([last_t_end;this_t]/1e6,[prop_output(i-1).population(end,:,j);prop_output(i).population(:,:,j)]*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
        end
        hold on;
    end
    hold off;
    %xlabel('Time (\mus)')
    ylabel('N_1 (%)');
    ylim([min(prop_output(1).population(:,:,j),[],1),max(prop_output(2).population(:,:,j),[],1)]*100);
    set(gca,'YColor','b');
    set(gca,'fontsize',20);
    if j == 1
        pos = get(gcf,'Position');
        set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
    end
    xlim([0,1/rep_rate*1e6]);
    %print(gcf,'all_pulses.pdf','-dpdf','-bestfit','-vector');

    % Add an inset for pulse's temporal profile
    ifig = figure;
    I_pulse = abs(prop_output(2).fields.forward(:,:,j)).^2; I_pulse = I_pulse/max(I_pulse);
    plot(t/1e3,I_pulse,'linewidth',2,'Color',[0.8510,0.3255,0.0980]);
    set(gca,'fontsize',12);
    xlabel('Time (ns)');
    ylabel('Power (norm.)');

    inset_fig = findobj(ifig,'Type','axes');
    h_inset = copyobj(inset_fig,fig);
    set(fig,'Units','normalized'); ax = get(fig,'Position');
    set(h_inset,'Position', [0.62 0.33 0.5*ax(4) 0.7*ax(4)])
    close(ifig);
    
    % Add an inset for population
    ifig = figure;
    plot(t/1e3,prop_output(2).population(:,:,j)*100,'linewidth',2,'Color','b');
    set(gca,'fontsize',12,'YColor','b','YAxisLocation','Right','Color','None');
    xlabel('');
    ylim([min(prop_output(2).population(:,:,j)),max(prop_output(2).population(:,:,j))]*100);
    ylabel('N_1 (%)');

    inset_fig2 = findobj(ifig,'Type','axes');
    h_inset2 = copyobj(inset_fig2,fig);
    set(fig,'Units','normalized'); ax = get(fig,'Position');
    set(h_inset2,'Position', [0.62 0.33 0.5*ax(4) 0.7*ax(4)])
    close(ifig);

    % Add fiber position
    iffig = figure;
    plot([0;prop_output(1).z(end)],[0.5;0.5],'linewidth',4,'Color','k');
    hold on;
    plot(prop_output(1).z(j),0.5,'Marker','o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','None');
    hold off;
    set(gca,'Box','Off','XTick',[],'YTick',[],'XColor','None','YColor','None','Color','None');

    fiber_fig = findobj(iffig,'Type','axes');
    h_fiber = copyobj(fiber_fig,fig);
    set(fig,'Units','normalized'); ax = get(fig,'Position');
    set(h_fiber,'Position', [0.2 0.65 0.5*ax(4) 0.7*ax(4)])
    close(iffig);
    text(0.15,0.88,'gain fiber','FontSize',20,'Units','Normalized');
    text(0.07,0.74,'in','FontSize',20,'Units','Normalized');
    text(0.33,0.74,'out','FontSize',20,'Units','Normalized');
    text(0.4,0.88,['\Delta N_1=',sprintf('%5.4f',abs(diff([min(prop_output(1).population(:,:,j),[],1),max(prop_output(1).population(:,:,j),[],1)]*100))),'%'],'FontSize',20,'Units','Normalized','Color','b');

    Frame(j) = getframe(fig);

    % Delete the inset for the next inset to plot on this figure
    delete(h_inset);
    delete(h_inset2);
    delete(h_fiber);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('Transient-gain');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);