%% Make Plots for ML-UAV paper, TRO journal: V, r, Zu,

close all
Environment = '1';
%Environment = '2';
%Environment = '3';

l_b = 0.8;
h_b = 0.25;
V_b = l_b*h_b^2;

t0=0;
t1=80-t0;
dt=1/500; %2ms
ind0=round(t0/dt)+1;
ind1=round(t1/dt);
%% Load and Prepare Data: 
%load('U');
%load('P_gyro_Noisy');
file_name = ['Monte Carlo\CM_1',Environment,'\SVCS\']; % new
load_parameters()
t=alpha.time(ind0:ind1,1);

V_l1_SVC =Xdot_b.data(ind0:ind1,1);
V_bar_l1_SVC =V_bar.data(ind0:ind1,1);
% V_w_l1_SVC =V_w.data(ind0:ind1,1);

X_b_l1_SVC = X_b.data(ind0:ind1,2);
zeta_w_l1_SVC = zeta_w.data(ind0:ind1,1);

r_l1_SVC = r.data(ind0:ind1,1);
r_bar_l1_SVC = r_bar.data(ind0:ind1,1);

z_u_l1_SVC = X_q.data(ind0:ind1,2)-10;
z_u_bar_l1_SVC = z_u_bar.data(ind0:ind1,1)-10;

active_l1_SVC =active.data(ind0:ind1,1);
coupled_est_l1_SVC =coupled_est.data(ind0:ind1,1);

alpha_l1_SVC = alpha.data(ind0:ind1,1);
alpha_bar_l1_SVC = alpha_bar.data(ind0:ind1,1);
switch_side_l1_SVC =switch_side.data(ind0:ind1,1);

theta_u_l1_SVC = 180/pi*theta_q.data(ind0:ind1,1);
theta_u_c_l1_SVC = 180/pi*theta_c.data(ind0:ind1,1);

T_l1_SVC =T.data(ind0:ind1,1);
T_bar_l1_SVC =T_bar.data(ind0:ind1,1);

Vol_im_l1_SVC =V_im.data(ind0:ind1,1);
Vol_im_bar_l1_SVC =V_im_bar.data(ind0:ind1,1);
Vol_im_bar_l1_SVC(1:5) = Vol_im_bar_l1_SVC(5,1)*ones(5,1);

u_1_l1_SVC = u_1c.data(ind0:ind1,1);
u_2_l1_SVC = u_2c.data(ind0:ind1,1);

% -------------------------------------------------------------------------
file_name = ['Monte Carlo\CM_2',Environment,'\SVCS\'];
load_parameters()

V_l2_SVC = Xdot_b.data(ind0:ind1,1);
V_bar_l2_SVC=V_bar.data(ind0:ind1,1);
% V_w_l2_SVC=V_w.data(ind0:ind1,1);

X_b_l2_SVC = X_b.data(ind0:ind1,2);
zeta_w_l2_SVC = zeta_w.data(ind0:ind1,1);

r_l2_SVC = r.data(ind0:ind1,1);
r_bar_l2_SVC = r_bar.data(ind0:ind1,1);

z_u_l2_SVC = X_q.data(ind0:ind1,2)-10;
z_u_bar_l2_SVC = z_u_bar.data(ind0:ind1,1)-10;

active_l2_SVC =active.data(ind0:ind1,1);
coupled_est_l2_SVC =coupled_est.data(ind0:ind1,1);

alpha_l2_SVC = alpha.data(ind0:ind1,1);
alpha_bar_l2_SVC = alpha_bar.data(ind0:ind1,1);
switch_side_l2_SVC =switch_side.data(ind0:ind1,1);

theta_u_l2_SVC = 180/pi*theta_q.data(ind0:ind1,1);
theta_u_c_l2_SVC = 180/pi*theta_c.data(ind0:ind1,1);

T_l2_SVC=T.data(ind0:ind1,1);
T_bar_l2_SVC=T_bar.data(ind0:ind1,1);

Vol_im_l2_SVC = V_im.data(ind0:ind1,1);
Vol_im_bar_l2_SVC = V_im_bar.data(ind0:ind1,1);
Vol_im_bar_l2_SVC(1:5) = Vol_im_bar_l2_SVC(5,1)*ones(5,1);

u_1_l2_SVC = u_1c.data(ind0:ind1,1);
u_2_l2_SVC = u_2c.data(ind0:ind1,1);
% -------------------------------------------------------------------------
% overlay_PID_results_l1_l2_compact()

time0=0:dt:(length(t)-1)*dt;

Dticks_time = 20;
ticks_time_end = t1-t0;
samp_red = 150; % plot at lower sample rate to reduce size

% -----------------------------------
%e_V_l1_PID = sum(abs(V_l1_PID - V_bar_l1_SVC))/(ind1-ind0)
e_V_l1_SVC = sum(abs(V_l1_SVC - V_bar_l1_SVC))/(ind1-ind0)
%e_V_l2_PID = sum(abs(V_l2_PID - V_bar_l2_SVC))/(ind1-ind0)
e_V_l2_SVC = sum(abs(V_l2_SVC - V_bar_l2_SVC))/(ind1-ind0)

ind_t10=round(10/dt);
ind_t43=round(43/dt);
ind_t57=round(57/dt);

V_bar_l2_SVC_active = [V_bar_l2_SVC(ind_t10:ind_t43);V_bar_l2_SVC(ind_t57:ind1)];
V_bar_l1_SVC_active = [V_bar_l1_SVC(ind_t10:ind_t43);V_bar_l1_SVC(ind_t57:ind1)];
%V_l2_PID_active = [V_l2_PID(ind_t10:ind_t43);V_l2_PID(ind_t57:ind1)];
V_l2_SVC_active = [V_l2_SVC(ind_t10:ind_t43);V_l2_SVC(ind_t57:ind1)];
%V_l1_PID_active = [V_l1_PID(ind_t10:ind_t43);V_l1_PID(ind_t57:ind1)];
V_l1_SVC_active = [V_l1_SVC(ind_t10:ind_t43);V_l1_SVC(ind_t57:ind1)];

%e_V_l2_PID_active = sum(abs(V_l2_PID_active - V_bar_l2_SVC_active))/length(V_bar_l2_SVC_active)
e_V_l2_SVC_active = sum(abs(V_l2_SVC_active - V_bar_l2_SVC_active))/length(V_bar_l2_SVC_active)
%e_V_l1_PID_active = sum(abs(V_l1_PID_active - V_bar_l1_SVC_active))/length(V_bar_l2_SVC_active)
e_V_l1_SVC_active = sum(abs(V_l1_SVC_active - V_bar_l1_SVC_active))/length(V_bar_l2_SVC_active)



%max(abs(V_l2_PID - V_bar_l1_SVC))
max(abs(V_l2_SVC - V_bar_l2_SVC))
%max(abs(V_l1_PID - V_bar_l1_SVC))
max(abs(V_l1_SVC - V_bar_l1_SVC))

%e_z_u_l2_PID = sum(abs(z_u_l2_PID - z_u_bar_l2_PID))/(ind1-ind0)
e_z_u_l2_SVC = sum(abs(z_u_l2_SVC - z_u_bar_l2_SVC))/(ind1-ind0)
%e_z_u_l1_PID = sum(abs(z_u_l1_PID - z_u_bar_l1_PID))/(ind1-ind0)
e_z_u_l1_SVC = sum(abs(z_u_l1_SVC - z_u_bar_l1_SVC))/(ind1-ind0)

%max(abs(z_u_l2_PID - z_u_bar_l2_PID))
max(abs(z_u_l2_SVC - z_u_bar_l2_SVC))
%max(abs(z_u_l1_PID - z_u_bar_l1_PID))
max(abs(z_u_l1_SVC - z_u_bar_l1_SVC))

% pause()


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% ------------------------- Plot V ----------------------------
close all
x0=9; 
y0=0.8; %position of the lower left side of the figure
width=3.85;
height=4.0;
figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

% -----------------------------Sub 1.1---------------------------------------
% sp_h = 0.1;
% sp_w = 0.35;
% sp_y_s = 0.05;
% sp_x_s = 0.05;
% sp_x0_c1 = 0.12;
% sp_x0_c2 = 0.05+sp_w;

% pos = [sp_x0_c1 0.86 sp_w sp_h];
ax01 = subplot(5,2,1);
% ax01 = subplot('Position',pos);

%%
hold on
%plot(time0(1:samp_red:end),V_l2_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),V_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),V_bar_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)

%plot(time0(1:samp_red:end),5*active_l1_SVC(1:samp_red:end),'b','LineWidth',0.8)

ylabel('$V$ (m/s)','Interpreter','latex')
% xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
    %'                                                                      ']})
axis([0 t1-t0 -5 5.8])

 % % 'XTick',0:%Dticks_time%:ticks_time_end,...
set(gca,...
'Units','normalized',...
'YTick',-5:2.5:10,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

% lgd=legend('PID','FSV','ref');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')


%
%
%
% 
%% ------------------------- Plot z_b ----------------------------

ax1 = subplot(5,2,2);

hold on
plot(time0(1:samp_red:end),X_b_l1_SVC(1:samp_red:end)-10,'r','LineWidth',0.8)
plot(time0(1:samp_red:end),zeta_w_l1_SVC(1:samp_red:end)-10,'k--','LineWidth',0.8)
%plot(time0(1:samp_red:end),X_b_l1_PID(1:samp_red:end)-10,'b','LineWidth',0.8)
%plot(time0(1:samp_red:end),zeta_w_l1_PID(1:samp_red:end)-10,'k:','LineWidth',0.8)

ylabel('$z$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(c)']})

axis([0 t1-t0 -2.5 2.5])

set(gca,...
'Units','normalized',...
'YTick',-3:1:3,...
'XTick',0:Dticks_time:ticks_time_end,...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

box on
grid on
set(gca,'gridlinestyle','-')

%% ------------------------- Plot r ----------------------------
% -----------------------------Sub 2.1---------------------------------------

ax03 = subplot(5,2,3);
% pos = [sp_x0_c1 0.75 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot(time0(1:samp_red:end),r_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),r_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),r_bar_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)

%plot(time0(1:samp_red:end),5*coupled_est_l1_SVC(1:samp_red:end),'b','LineWidth',0.8)

ylabel('$r$ (m)','Interpreter','latex')
%xlabel({'Time (s)';...
%    ['(c) case 1: current only.']})%...
    %'                                                                      ']})
axis([0 t1-t0 3 8])

set(gca,...
'Units','normalized',...
'YTick',1:3:20,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')


% lgd=legend('PID','FSV','ref');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')
% ------------------------------------------------
%% ----------------------------- Plot z_u   -------------------------------
ax05 = subplot(5,2,5);
% pos = [sp_x0_c1 0.65 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
%plot([1,70],[5 5],'--k','LineWidth',0.8)
%plot(time0(1:samp_red:end),z_u_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),z_u_bar_l1_SVC(1:samp_red:end),'k-.','LineWidth',0.8)
%plot(time0(1:samp_red:end),z_u_bar_l1_PID(1:samp_red:end),'k:','LineWidth',0.8)

ylabel('$z_{\mathrm{u}}$ (m)','Interpreter','latex')
%xlabel('Time (s)')
%xlabel({'Time (s)';...
%    ['(e) case 1: current only.']})

axis([0 t1-t0 2 8])

set(gca,...
'Units','normalized',...
'YTick',3:2:10,...
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

% ----------------------------------------------------------

% #########################################################################

% #########################################################################

% ------------------------------------------------

%% ------------------------- Plot alpha theta ----------------------------

% -----------------------------Sub 1.1---------------------------------------
ax04 = subplot(5,2,4);
% pos = [sp_x0_c2 0.18 sp_w sp_h];
% ax01 = subplot('Position',pos);

hold on
plot(time0(1:samp_red:end),alpha_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
%plot(time0(1:samp_red:end),alpha_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)
plot(time0(1:samp_red:end),alpha_bar_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)
%ss_l2 = plot(time0(1:samp_red:end),90*switch_side_l1_SVC(1:samp_red:end),'k','LineWidth',0.8);

ylabel('\alpha (\circ)')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 25 160])

set(gca,...
'Units','normalized',...
'YTick',0:45:180,...
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

% -----------------------------Sub 2.1---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 1:

% -----------------------------Sub 2.1---Theta-------------------------------
ax06 = subplot(5,2,6);
box on
grid on
set(gca,'gridlinestyle','-')

hold on
%plot(time0(1:samp_red:end),theta_u_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),theta_u_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),theta_u_c_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)
% plot(time0(1:samp_red:end),theta_u_c_l1_PID(1:samp_red:end),'k-.','LineWidth',0.8)

ylabel('\theta_{u} (\circ)')
%xlabel({'Time (s)';...
%    ['(c) case 1: current only.']})%,...
    %'                                       ']})

axis([0 t1-t0 -50 50])

set(gca,...
'Units','normalized',...
'YTick',-40:40:80,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------- Plot T and Vol ----------------------------

% -----------------------------Sub 1.1---------------------------------------

ax07 = subplot(5,2,7);

hold on
%plot(time0(1:samp_red:end),T_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),T_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)
plot(time0(1:samp_red:end),T_bar_l1_SVC(1:samp_red:end),'k--','LineWidth',0.8)

ylabel('$T$ (N)','Interpreter','latex')
%xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
axis([0 t1-t0 -5 100])

set(gca,...
'Units','normalized',...
'YTick',0:25:100,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

box on
grid on
set(gca,'gridlinestyle','-')

% ------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%------------------------- Plot  Vol ----------------------------

% -------------------------- 2.1
ax08 = subplot(5,2,8);
hold on

%plot(time0(1:samp_red:end),Vol_im_l1_PID(1:samp_red:end)/V_b,'b','LineWidth',0.8)
plot(time0(1:samp_red:end),Vol_im_l1_SVC(1:samp_red:end)/V_b,'r','LineWidth',0.8)
% plot(time0(1:samp_red:end),Vol_im_bar_l1_SVC(1:samp_red:end)/V_b,'--k','LineWidth',0.8)

ylabel('$V_{\mathrm{im}}/V_{\mathrm{b}}$','Interpreter','latex') %curlyvee_
%xlabel({'Time (s)';...
%    ['(c) case 1: current only.']})%...
    %'                                            ']})

axis([0 t1-t0 0.0 0.5])

set(gca,...
'Units','normalized',...
'YTick',0.0:0.25:1,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
xticklabels('')

% lgd=legend('PC estimated mean','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')

% #########################################################################


%% ------------------------- Plot u1 u2 ----------------------------


% -----------------------------Sub 1.1---------------------------------------

ax09 = subplot(5,2,9);

hold on
%plot(time0(1:samp_red:end),u_1_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),u_1_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)

ylabel('u_{1c} (N)')
xlabel({'Time (s)';...
    ['                                                      ',...
    '(a) S1: mini-size multirotor with cable length of 7 m.']})
axis([0 t1-t0 0 126])

set(gca,...
'Units','normalized',...
'YTick',0:25:200,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')

% lgd=legend('PC ref','PID','PC');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')

% -----------------------------Sub 2.1---------------------------------------
ax10 = subplot(5,2,10);

hold on
%plot(time0(1:samp_red:end),u_2_l1_PID(1:samp_red:end),'b','LineWidth',0.8)
plot(time0(1:samp_red:end),u_2_l1_SVC(1:samp_red:end),'r','LineWidth',0.8)

ylabel('u_{2c} (N.m)')
xlabel('Time (s)')

axis([0 t1-t0 -3.5 3.5])

set(gca,...
'Units','normalized',...
'YTick',-3:3:3,...
'XTick',0:Dticks_time:ticks_time_end,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')

% lgd=legend('desired','actual');
% set(lgd,'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',7,...
% 'FontName','Times',...
% 'Location','SouthEast')

box on
grid on
set(gca,'gridlinestyle','-')


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%      -----------------------------
% print -depsc2 XYZ_ICS_Sim.eps
% save as pdf:
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

%print -dpdf -painters [Plots/all_states_S1_E2_0]
print(['Plots/all_states_S1_E',Environment,'_0'],'-dpdf')
print(['Plots/all_states_S1_E',Environment,'_0'],'-dpng','-r500')

% -----------------------------------------------

plot_S2_Ej;
% -----------------------------------------------------------------------
%      -----------------------------
% print -depsc2 XYZ_ICS_Sim.eps
% save as pdf:
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters Plots/all_states_S2_E2_0
print(['Plots/Monte Carlo/all_states_S2_E',Environment,'_0'],'-dpdf')
print(['Plots/Monte Carlo/all_states_S2_E',Environment,'_0'],'-dpng','-r500')

