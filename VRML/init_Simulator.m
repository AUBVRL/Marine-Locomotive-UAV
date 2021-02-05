%% Initialize simulator:

T = 1/30; % Simulation time step [s]
Lc_0 = 7; % cable length [m]
g = 9.81; % gravity constant [m/s^2]

H = 10;
alpha_bar_0 = pi/4;     % cable angle
Delta_h = 0;
X0_b = [0,H-Delta_h];
zu_bar = X0_b(2) + Lc_0*sin(alpha_bar_0);


%% Wave characteristics
% AWTEC 2016: beirut, H = 2*A = 0.3 m, T = 3.5 s
% Wiki: fully developped, H =0.27 m, T = 3 s
% Wiki: fully developped, H =1.5 m, T = 5.7 s
% Wiki: fully developped, H =4.1 m, T = 8.6 s
% Wiki: Interpolate: H =2 m, T = 6.3 s
nu_w = 1.787*10^-6;   % kinematic viscosity of water (m2/s)
A_w = 1*1.5/2;          % wave height (m)
w_dir = 1;           % wave direction (1 for 0 deg, -1 for 180 deg)
H = 10;               % mean water level (m)
T_w = 5.7;            % wave period (s)
omega_w = 2*pi/T_w;   % wave circular frequency (rad/s)
k_w = omega_w^2/g;    % wave number (rad/m)
epsilon_w = 0;        % random phase angle (rad)
C_w = omega_w/k_w;    % wave speed (m/s)
lambda_w = C_w*T_w;   % wave length (m)

%
A_w2 = 1*0.27/2;        % wave height (m)
w_dir2 = 1;           % wave direction (1 for 0 deg, -1 for 180 deg)
T_w2 = 3;            % wave period (s)
omega_w2 = 2*pi/T_w2;   % wave circular frequency (rad/s)
k_w2 = omega_w2^2/g;    % wave number (rad/m)
epsilon_w2 = pi;        % random phase angle (rad)
%