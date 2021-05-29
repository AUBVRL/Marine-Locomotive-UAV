% Monte carlo
rho_w = 1000;
g = 9.81;

nu_w = 1.787*10^-6;   % kinematic viscosity of water (m2/s)
% H = 10;               % mean water level (m)
% T_w = 5.7;            % wave period (s)
% omega_w = 2*pi/T_w;   % wave circular frequency (rad/s)
% k_w = omega_w^2/g;    % wave number (rad/m)
% epsilon_w = 0;        % random phase angle (rad)
% C_w = omega_w/k_w;    % wave speed (m/s)
% lambda_w = C_w*T_w;   % wave length (m)

h_b = 0.25;            % buoy side (m), cuboid 
A_b = h_b^2;           % (m^2)
l_b = 0.8;             % buoy length (m)
Vol_b = A_b*l_b;       % buoy volume (m)

% for the buoy to be one quarter immersed:
m_b = rho_w*Vol_b/4;   % buoy mass (kg)
J_b = m_b*h_b*l_b/12; 

a_11 = 0.05*m_b;        % added mass
a_33 = m_b;
b_11 = 0;               % added dampind

C_S = [5;9]*10^-3;       % skin friction constant
D_s = [4;0.3*3];         % sfkin friction coefficient (viscous)
%D_x = b_11 + D_s(1);     % total drag coefficient


A_im = (m_b)/rho_w/l_b; % immersed frontal area (m^2) (different from whetted area)
Delta_h = (A_im/A_b-0.5)*h_b; % z_b-zeta, buoy center above water surface 

%% Frequency analysys:
A_c = l_b*h_b;
D_z0 = 55;

%Gb_s = tf(omega_b^2,[1 2*zeta_b*omega_b omega_b^2])
omega_b = sqrt((rho_w*g*A_c/(m_b+a_33)));


%A_wave = [0.27;0.43;0.61;0.85;1.2;1.5;1.9;3.3]/2;
%T_w = [3;3.4;4;4.6;5;5.7;6.3;8];            % wave period (s)
A_wave = [0.27;0.61;1.2;1.5;3.3]/2;
T_w = [3;4;5;5.7;8];            % wave period (s)

omega_w = 2*pi/T_w;   % wave circular frequency (rad/s)
V_bar = [-20:0.1:20]';

Vol_im_0 = m_b/rho_w;

for i = 1:length(A_wave)
    T_i = T_w(i);            % wave period (s)
    A_i = A_wave(i);
    omega_w = 2*pi/T_i;   % wave circular frequency (rad/s)
    k_w = omega_w^2/g;    % wave number (rad/m)
    w_dir = 1;           % wave direction (1 for 0 deg, -1 for 180 deg)
    
    u_s(i) = w_dir*A_i^2*omega_w*k_w*exp(2*k_w*(-Delta_h/2)); % stocks drift (m/s)

    omega_h = omega_w;
    b_33 = 2*m_b*omega_h;   
    D_z = b_33 + D_s(2);
    zeta_b = D_z/(2*sqrt((m_b+a_33)*rho_w*g*A_c));
        
            
    
    for j = 1:length(V_bar)
        alpha_bar = pi/4*sign(V_bar(j));     % cable angle
        V_bar_j = V_bar(j) - u_s(i);

        Re = abs(V_bar_j*l_b/nu_w); % Reynolds number
        %Fn(i,1) = V_b/sqrt(g*l_b); % froud number, describes the velocity 
        %K_cor = 4.44/(515.38*0.3048^2*0.3048^2); % unit conversion correction,
        %equal to unity
        C_S = 0.075/(log10(Re)-2)^2; % skin friction coeff
        % C_S(i,1) = 0.003;
        % F_D(i,1) = C_F*(0.5*rho_w*V_b^2)*(2*R_b*theta*l_b); % cylinder
        D_s(i,j) = C_S*(0.5*rho_w*abs(V_bar_j))*(l_b*h_b + l_b*h_b/4*0.0 + l_b*h_b/4*0.5); % Cuboid    
        D_x = b_11 + D_s(i,j);     % total drag coefficient


        Vol_im(i,j) = m_b/rho_w - (V_bar_j/rho_w/g)*(D_x*tan(alpha_bar));
        dh_bar(i,j) = Vol_im(i,j)/A_c;

        omega_e = abs(omega_w - w_dir*omega_w^2*V_bar_j/g); % -0.508 1.72
        OMEGA = omega_e/omega_b;
        dh_amp(i,j) = A_i*(1/sqrt((1 - OMEGA^2)^2+(2*zeta_b*OMEGA)^2)-1); % cm
        dVol_amp(i,j) = dh_amp(i,j)*A_c;
        Vol_im_ratio(i,j) = dVol_amp(i,j)/Vol_im(i,j);
    end
end

%% Plot
close all
x0=9; 
y0=0.8; %position of the lower left side of the figure
width=3.85;
height=2.0;
figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

hold on
plot(V_bar,dh_bar(1,:),'k','LineWidth',1.2)
for k = 1:length(A_wave)
    plot(V_bar,dh_amp(k,:))
end
grid on


ylabel('$\Delta h_{\mathrm{amp},n}$ (m)','Interpreter','latex')
% xlabel({'Time (s)';...
%    ['(a) case 1: current only.']})%...
    %'                                                                      ']})
xlabel('$V$ (m/s)','Interpreter','latex')
title('$S_{\bar{V},n}^{\mathrm{fo}}$','Interpreter','latex')
axis([-20 20 0 0.07])

 % % 'XTick',0:%Dticks_time%:ticks_time_end,...
set(gca,...
'Units','normalized',...
'YTick',0:0.01:0.07,...
'XTick',-20:5:20,... 
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',8,...
'FontName','Times')
%xticklabels('')

lgd=legend('$\bar{h}_{\mathrm{im}}$','n=1','n=2','n=3','n=4','n=5');
set(lgd,'FontUnits','points',...
'interpreter','latex',...
'FontSize',8,...
'FontName','Times',...
'Location','Best')

box on
grid on
set(gca,'gridlinestyle','-')

%% 
% figure(2)
% hold on
% plot(V_bar,dh_bar(2,:))
% plot(V_bar,dh_amp(2,:))
% grid on
% 
% figure(6)
% hold on
% plot(V_bar,dh_bar(6,:))
% plot(V_bar,dh_amp(6,:))
% grid on
% 
% figure(7)
% hold on
% plot(V_bar,dh_bar(7,:))
% plot(V_bar,dh_amp(7,:))
% grid on
 
% -----------------------------------------------------------------------
%      -----------------------------
% print -depsc2 XYZ_ICS_Sim.eps
% save as pdf:
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters Plots/dynamic_amplification