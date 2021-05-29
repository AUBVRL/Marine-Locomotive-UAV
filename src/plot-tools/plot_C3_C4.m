close all

%% Load and Prepare Data: 
l_b = 0.8;
h_b = 0.25;
V_b = l_b*h_b^2;

%%
t0=20;
t1=100-t0;
dt=1/500; %2ms
ind0_C3=round(t0/dt)+1;
ind1_C3=round(t1/dt);
H=10;

file_name = 'C3\';
load_parameters()
V_C3 =Xdot_b.data(ind0_C3:ind1_C3,1);
V_bar_C3 =V_bar.data(ind0_C3:ind1_C3,1);
T_C3 =T.data(ind0_C3:ind1_C3,1);
Vol_im_C3 =V_im.data(ind0_C3:ind1_C3,1);
z_u_C3 = X_q.data(ind0_C3:ind1_C3,2)-10;
z_u_bar_C3 = z_u_bar.data(ind0_C3:ind1_C3,1)-10;
alpha_C3 = alpha.data(ind0_C3:ind1_C3,1);
alpha_bar_C3 = alpha_bar.data(ind0_C3:ind1_C3,1);

X_b_C3 = X_b.data(ind0_C3:ind1_C3,2);
zeta_w_C3 = zeta_w.data(ind0_C3:ind1_C3,1);

t_C3=Xdot_b.time(ind0_C3:ind1_C3,1);
time0=0:dt:(length(t_C3)-1)*dt;

%%
t0=20;
t1=100-t0;
dt=1/500; %2ms
ind0_C4=round(t0/dt)+1;
ind1_C4=round(t1/dt);

file_name = 'C4\A165cm\';
load_parameters()
V_C4 =Xdot_b.data(ind0_C4:ind1_C4,1);
V_bar_C4 =V_bar.data(ind0_C4:ind1_C4,1);
T_C4 =T.data(ind0_C4:ind1_C4,1);
Vol_im_C4 =V_im.data(ind0_C4:ind1_C4,1);
z_u_C4 = X_q.data(ind0_C4:ind1_C4,2)-10;
z_u_bar_C4 = z_u_bar.data(ind0_C4:ind1_C4,1)-10;
alpha_C4 = alpha.data(ind0_C4:ind1_C4,1);
alpha_bar_C4 = alpha_bar.data(ind0_C4:ind1_C4,1);
X_b_C4 = X_b.data(ind0_C4:ind1_C4,2);
zeta_w_C4 = zeta_w.data(ind0_C4:ind1_C4,1);

t_C4=Xdot_b.time(ind0_C4:ind1_C4,1);
time0=0:dt:(length(t_C4)-1)*dt;

%%
Dticks_time = 20;
ticks_time_end = t1-t0;
ticks_time_end = t1-t0;

samp_red = 50; % plot at lower sample rate to reduce size

% -----------------------------------
%% ------------------------- Plot alpha theta V----------------------------

x0=9; 
y0=4; %position of the lower left side of the figure
width=3.85;
height=4.5;
figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

% -----------------------------------
% -----------------------------------
% -----------------------------------

ax1 = subplot(6,2,1);

hold on
plot(time0(1:samp_red:end),V_C3(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),V_bar_C3(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$V$ (m/s)','Interpreter','latex')
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(a)']})

axis([0 t1-t0 0 15])

set(gca,...
'Units','normalized',...
'YTick',0:5:20,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

% lgd=legend('actual','desired');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')

% lgd=legend('wave-free','forward-waves','backward-waves');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

% -----------------------------------

% ------------------------------------------------------------
% ----------------------------------------------------------

ax1 = subplot(6,2,2);

hold on
plot(time0(1:samp_red:end),Vol_im_C3(1:samp_red:end)/V_b,'r','LineWidth',0.8)


ylabel('$V_{\mathrm{im}}/V_{\mathrm{b}}$','Interpreter','latex') %curlyvee_
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(e)']})

axis([0 t1-t0 -0.1 1.1])

set(gca,...
'Units','normalized',...
'YTick',0:0.5:10,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

% lgd=legend('ref','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')


%% Plot Zu
ax05 = subplot(6,2,3);
% pos = [sp_x0_c1 0.65 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot([1,70],[5 5],'--k','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_C3(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_bar_C3(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$z_{\mathrm{u}}$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(e) case 1: current only.']})

axis([0 t1-t0 4.8 7.1])

set(gca,...
'Units','normalized',...
'YTick',1:1:10,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

% lgd=legend('actual','desired');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')


% lgd=legend('wave-free','forward-waves','backward-waves');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

% ----------------------------------------------------------
%% ------------------------- Plot alpha ----------------------------

% -----------------------------Sub 1.1---------------------------------------
ax04 = subplot(6,2,4);
% pos = [sp_x0_c2 0.18 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
plot(time0(1:samp_red:end),alpha_C3(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_C3(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('\alpha (\circ)')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 30 60])

set(gca,...
'Units','normalized',...
'YTick',30:10:60,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

% lgd=legend('PC ref','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')

% ------------------------------------------------
% ------------------------------------------------
% 
ax1 = subplot(6,2,5:6);

hold on
plot(time0(1:samp_red:end),X_b_C3(1:samp_red:end)-10,'r','LineWidth',0.8)
plot(time0(1:samp_red:end),zeta_w_C3(1:samp_red:end)-10,'c','LineWidth',0.6)

ylabel('$z$ (m)','Interpreter','latex')
%xlabel('Time (s)')
xlabel({'Time (s)';...
    ['(a) case 3: high frequency, head seas.']})

axis([0 t1-t0 -0.2 1.5])

set(gca,...
'Units','normalized',...
'YTick',0:0.5:2,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')

box on
grid on
set(gca,'gridlinestyle','-')

% ------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------



% -----------------------------------
% -----------------------------------
ax1 = subplot(6,2,7);

hold on
plot(time0(1:samp_red:end),V_C4(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),V_bar_C4(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$V$ (m/s)','Interpreter','latex')
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(b)']})

axis([0 t1-t0 0 15])

set(gca,...
'Units','normalized',...
'YTick',0:5:15,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

box on
grid on
set(gca,'gridlinestyle','-')

% -------------------------------------------
% -------------------------------------------
% -------------------------------------------
% -------------------------------------------

ax1 = subplot(6,2,8);

hold on
plot(time0(1:samp_red:end),Vol_im_C4(1:samp_red:end)/V_b,'r','LineWidth',0.8)


ylabel('$V_{\mathrm{im}}/V_{\mathrm{b}}$','Interpreter','latex') %curlyvee_
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(f)']})

axis([0 t1-t0 -0.1 1.1])

set(gca,...
'Units','normalized',...
'YTick',0:0.5:10,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

% lgd=legend('ref','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')

% ------------------------------------------------

%% Plot Zu
ax05 = subplot(6,2,9);
% pos = [sp_x0_c1 0.65 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot([1,70],[5 5],'--k','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_C4(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_bar_C4(1:samp_red:end),'k-.','LineWidth',0.8)

ylabel('$z_{\mathrm{u}}$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(e) case 1: current only.']})

axis([0 t1-t0 4.5 9.5])

set(gca,...
'Units','normalized',...
'YTick',5:2:10,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

% lgd=legend('actual','desired');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')


% lgd=legend('wave-free','forward-waves','backward-waves');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

% ----------------------------------------------------------
%% ------------------------- Plot alpha ----------------------------

ax04 = subplot(6,2,10);
% pos = [sp_x0_c2 0.18 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
plot(time0(1:samp_red:end),alpha_C4(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_C4(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('\alpha (\circ)')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 10 82])

set(gca,...
'Units','normalized',...
'YTick',25:25:100,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

% lgd=legend('PC ref','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')

% ------------------------------------------------
%% Zeta z_b

ax1 = subplot(6,2,11:12);

hold on
plot(time0(1:samp_red:end),X_b_C4(1:samp_red:end)-10,'r','LineWidth',0.8)
plot(time0(1:samp_red:end),zeta_w_C4(1:samp_red:end)-10,'c','LineWidth',0.6)

ylabel('$z$ (m)','Interpreter','latex')
%xlabel('Time (s)')
xlabel({'Time (s)';...
    ['(b) case 4: high amplitude, head seas.']})

axis([0 t1-t0 -2.5 3.5])

set(gca,...
'Units','normalized',...
'YTick',-1.5:1.5:4.5,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')

box on
grid on
set(gca,'gridlinestyle','-')
% ------------------------------------------------
% ------------------------------------------------------------

% ax1 = subplot(4,2,7);
% 
% hold on
% plot(time0(1:samp_red:end),T_C3(1:samp_red:end),'r','LineWidth',0.8)
% 
% ylabel('$T$ (N)','Interpreter','latex')
% %xlabel('Time (s)')
% xlabel({'Time (s)';...
%     '(g)';...
%     'Case 3: positive velocity.'})%...
% 
% axis([0 t1-t0 -5 150])
% 
% set(gca,...
% 'Units','normalized',...
% 'YTick',0:50:150,...
% 'XTick',0:Dticks_time:ticks_time_end,...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',8,...
% 'FontName','Times')
% 
% box on
% grid on
% set(gca,'gridlinestyle','-')

% ------------------------------------------------

% ------------------------------------------------
% ------------------------------------------------------------

% ax1 = subplot(4,2,8);
% 
% hold on
% plot(time0(1:samp_red:end),T_C4(1:samp_red:end),'r','LineWidth',0.8)
% 
% ylabel('$T$ (N)','Interpreter','latex')
% %xlabel('Time (s)')
% xlabel({'Time (s)';...
%     '(h)';...
%     'Case 4: negative velocity.'})%...
% 
% axis([0 t1-t0 -5 150])
% 
% set(gca,...
% 'Units','normalized',...
% 'YTick',0:50:150,...
% 'XTick',0:Dticks_time:ticks_time_end,...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',8,...
% 'FontName','Times')
% 
% box on
% grid on
% set(gca,'gridlinestyle','-')

% ------------------------------------------------

%      -----------------------------
% print -depsc2 XYZ_ICS_Sim.eps
% save as pdf:
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
 print -dpdf -painters Plots/C3_C4
