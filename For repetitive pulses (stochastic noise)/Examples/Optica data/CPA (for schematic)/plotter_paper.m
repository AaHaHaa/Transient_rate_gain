figure;
yyaxis left;
plot(t_woPulse,prop_output(1).population(:,:,end)*100,'linewidth',2,'Color','b');
ylabel('N_1 (%)');
yyaxis right;
plot(t_woPulse,abs(prop_output(1).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
print(gcf,'Aout_1_wASE.pdf','-dpdf','-bestfit');

% Just for schematic visualization purpose
% So I add more ASE noise with the signal to show that ASE is also in this
% window
A1 = interp1(abs(prop_output(1).fields.forward(:,:,end)).^2,linspace(1,1024,65536)');
A2 = abs(prop_output(2).fields.forward(:,:,end)).^2/1e3;
idx = A2 < max(A2)/10;
A2(idx) = A1(idx)*2 + A2(idx);
figure;
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2,'Color','b');
ylabel('N_1 (%)');
yyaxis right;
plot(t,A2,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
print(gcf,'Aout_2_wASE.pdf','-dpdf');