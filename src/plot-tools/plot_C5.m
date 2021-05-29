close all

%% Load and Prepare Data: 
l_b = 0.8;
h_b = 0.25;
V_b = l_b*h_b^2;

%%
t0=0;
t1=30-t0;
dt=1/500; %2ms
ind0_C5=round(t0/dt)+1;
ind1_C5=round(t1/dt);
ind1_C5_CBNC=round(8/dt);

H=10;

file_name = 'C5/drop/SVCS/';
load_parameters()
V_C5_drop_SVCS =Xdot_b.data(ind0_C5:ind1_C5,1);
V_bar_C5_drop_SVCS =V_bar.data(ind0_C5:ind1_C5,1);
T_C5_drop_SVCS =T.data(ind0_C5:ind1_C5,1);
Vol_im_C5_drop_SVCS =V_im.data(ind0_C5:ind1_C5,1);
z_u_C5_drop_SVCS = X_q.data(ind0_C5:ind1_C5,2)-10;
z_u_bar_C5_drop_SVCS = z_u_bar.data(ind0_C5:ind1_C5,1)-10;
alpha_C5_drop_SVCS = alpha.data(ind0_C5:ind1_C5,1);
alpha_bar_C5_drop_SVCS = alpha_bar.data(ind0_C5:ind1_C5,1);
X_b_C5_drop_SVCS = X_b.data(ind0_C5:ind1_C5,2);
zeta_w_C5_drop_SVCS = zeta_w.data(ind0_C5:ind1_C5,1);
r_C5_drop_SVCS = r.data(ind0_C5:ind1_C5,1);
r_bar_C5_drop_SVCS = r_bar.data(ind0_C5:ind1_C5,1);

t_C5=Xdot_b.time(ind0_C5:ind1_C5,1);
time0=0:dt:(length(t_C5)-1)*dt;

file_name = 'C5/drop/CBNC/';
load_parameters()
V_C5_drop_CBNC =Xdot_b.data(ind0_C5:ind1_C5,1);
V_bar_C5_drop_CBNC =V_bar.data(ind0_C5:ind1_C5,1);
T_C5_drop_CBNC =T.data(ind0_C5:ind1_C5,1);
Vol_im_C5_drop_CBNC =V_im.data(ind0_C5:ind1_C5,1);
z_u_C5_drop_CBNC = X_q.data(ind0_C5:ind1_C5,2)-10;
z_u_bar_C5_drop_CBNC = z_u_bar.data(ind0_C5:ind1_C5,1)-10;
alpha_C5_drop_CBNC = alpha.data(ind0_C5:ind1_C5,1);
alpha_bar_C5_drop_CBNC = alpha_bar.data(ind0_C5:ind1_C5,1);
r_C5_drop_CBNC = r.data(ind0_C5:ind1_C5,1);
r_bar_C5_drop_CBNC = r_bar.data(ind0_C5:ind1_C5,1);

X_b_C5_drop_CBNC = X_b.data(ind0_C5:ind1_C5,2);
zeta_w_C5_drop_CBNC = zeta_w.data(ind0_C5:ind1_C5,1);

%%
file_name = 'C5/moving/SVCS/';
load_parameters()
V_C5_moving_SVCS =Xdot_b.data(ind0_C5:ind1_C5,1);
V_bar_C5_moving_SVCS =V_bar.data(ind0_C5:ind1_C5,1);
T_C5_moving_SVCS =T.data(ind0_C5:ind1_C5,1);
Vol_im_C5_moving_SVCS =V_im.data(ind0_C5:ind1_C5,1);
z_u_C5_moving_SVCS = X_q.data(ind0_C5:ind1_C5,2)-10;
z_u_bar_C5_moving_SVCS = z_u_bar.data(ind0_C5:ind1_C5,1)-10;
alpha_C5_moving_SVCS = alpha.data(ind0_C5:ind1_C5,1);
alpha_bar_C5_moving_SVCS = alpha_bar.data(ind0_C5:ind1_C5,1);
X_b_C5_moving_SVCS = X_b.data(ind0_C5:ind1_C5,2);
zeta_w_C5_moving_SVCS = zeta_w.data(ind0_C5:ind1_C5,1);
r_C5_moving_SVCS = r.data(ind0_C5:ind1_C5,1);
r_bar_C5_moving_SVCS = r_bar.data(ind0_C5:ind1_C5,1);

file_name = 'C5/moving/CBNC/';
load_parameters()
V_C5_moving_CBNC =Xdot_b.data(ind0_C5:ind1_C5,1);
V_bar_C5_moving_CBNC =V_bar.data(ind0_C5:ind1_C5,1);
T_C5_moving_CBNC =T.data(ind0_C5:ind1_C5,1);
Vol_im_C5_moving_CBNC =V_im.data(ind0_C5:ind1_C5,1);
z_u_C5_moving_CBNC = X_q.data(ind0_C5:ind1_C5,2)-10;
z_u_bar_C5_moving_CBNC = z_u_bar.data(ind0_C5:ind1_C5,1)-10;
alpha_C5_moving_CBNC = alpha.data(ind0_C5:ind1_C5,1);
alpha_bar_C5_moving_CBNC = alpha_bar.data(ind0_C5:ind1_C5,1);
X_b_C5_moving_CBNC = X_b.data(ind0_C5:ind1_C5,2);
zeta_w_C5_moving_CBNC = zeta_w.data(ind0_C5:ind1_C5,1);
r_C5_moving_CBNC  = r.data(ind0_C5:ind1_C5,1);
r_bar_C5_moving_CBNC  = r_bar.data(ind0_C5:ind1_C5,1);

%%
Dticks_time = 10;
ticks_time_end = t1-t0;
ticks_time_end = t1-t0;

samp_red = 50; % plot at lower sample rate to reduce size

% -----------------------------------
%% ------------------------- Plot alpha theta V----------------------------

x0=9; 
y0=4; %position of the lower left side of the figure
width=3.85;
height=3.0;
figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

% -----------------------------------
% -----------------------------------
% -----------------------------------

ax1 = subplot(4,2,1);

hold on
plot(time0(1:samp_red:end),V_C5_drop_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:ind1_C5_CBNC),V_C5_drop_CBNC(1:samp_red:ind1_C5_CBNC),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),V_bar_C5_drop_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$V$ (m/s)','Interpreter','latex')
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(a)']})

axis([0 t1-t0 -6 5])

set(gca,...
'Units','normalized',...
'YTick',-5:2.5:5,...
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

ax1 = subplot(4,2,4);

hold on
%plot(time0(1:samp_red:end),r_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),r_C5_drop_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:ind1_C5_CBNC),r_C5_drop_CBNC(1:samp_red:ind1_C5_CBNC),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),r_bar_C5_drop_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

%plot(time0(1:samp_red:end),5*coupled_est_l1_SVC(1:samp_red:end),'b','LineWidth',0.8)

ylabel('$r$ (m)','Interpreter','latex')
%xlabel({'Time (s)';...
%    ['(c) case 1: current only.']})%...
    %'                                                                      ']})
axis([0 t1-t0 3 8])

set(gca,...
'Units','normalized',...
'YTick',2:2:20,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

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
ax05 = subplot(4,2,3);
% pos = [sp_x0_c1 0.65 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot([1,70],[5 5],'--k','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_C5_drop_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:ind1_C5_CBNC),z_u_C5_drop_CBNC(1:samp_red:ind1_C5_CBNC),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_bar_C5_drop_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$z_{\mathrm{u}}$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(e) case 1: current only.']})

axis([0 t1-t0 0 8])

set(gca,...
'Units','normalized',...
'YTick',0:2:10,...
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
ax04 = subplot(4,2,2);
% pos = [sp_x0_c2 0.18 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
plot(time0(1:samp_red:end),alpha_C5_drop_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:ind1_C5_CBNC),alpha_C5_drop_CBNC(1:samp_red:ind1_C5_CBNC),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_C5_drop_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('\alpha (\circ)')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 0 90])

set(gca,...
'Units','normalized',...
'YTick',0:30:90,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

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



% -----------------------------------
% -----------------------------------
ax1 = subplot(4,2,5);

hold on
plot(time0(1:samp_red:end),V_C5_moving_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),V_C5_moving_CBNC(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),V_bar_C5_moving_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$V$ (m/s)','Interpreter','latex')
%xlabel('Time (s)')
% xlabel({'Time (s)';...
%     ['(b)']})

axis([0 t1-t0 -6 5])

set(gca,...
'Units','normalized',...
'YTick',-10:2.5:5,...
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

ax1 = subplot(4,2,8);

hold on
%plot(time0(1:samp_red:end),r_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),r_C5_moving_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),r_C5_moving_CBNC(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),r_bar_C5_moving_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

%plot(time0(1:samp_red:end),5*coupled_est_l1_SVC(1:samp_red:end),'b','LineWidth',0.8)

ylabel('$r$ (m)','Interpreter','latex')
%xlabel({'Time (s)';...
%    ['(c) case 1: current only.']})%...
    %'                                                                      ']})
axis([0 t1-t0 3 8])

set(gca,...
'Units','normalized',...
'YTick',0:2:20,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

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
ax05 = subplot(4,2,7);
% pos = [sp_x0_c1 0.65 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot([1,70],[5 5],'--k','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_C5_moving_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_C5_moving_CBNC(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_bar_C5_moving_SVCS(1:samp_red:end),'k-.','LineWidth',0.8)

ylabel('$z_{\mathrm{u}}$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(e) case 1: current only.']})

axis([0 t1-t0 0 8])

set(gca,...
'Units','normalized',...
'YTick',0:2:10,...
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

ax04 = subplot(4,2,6);
% pos = [sp_x0_c2 0.18 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
plot(time0(1:samp_red:end),alpha_C5_moving_SVCS(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_C5_moving_CBNC(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_C5_moving_SVCS(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('\alpha (\circ)')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 10 95])

set(gca,...
'Units','normalized',...
'YTick',25:25:100,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

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

%      -----------------------------
% print -depsc2 XYZ_ICS_Sim.eps
% save as pdf:
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
 print -dpdf -painters Plots/C5_0
