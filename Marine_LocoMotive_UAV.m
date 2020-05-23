
T=1/500;              % sec, Simulation step time.
g = 9.81;             % gravity constant (m/s2)

%% Initalize:
alpha_bar = pi/4;     % cable angle
V_bar = 20;            % desired mean velocity
epsilon_1 = 6;        % cable tension lower bound (N)
%% Wave characteristics
% AWTEC 2016: beirut, H = 2*A = 0.3 m, T = 3.5 s
% Wiki: fully developped, H =0.27 m, T = 3 s
% Wiki: fully developped, H =1.5 m, T = 5.7 s
% Wiki: fully developped, H =4.1 m, T = 8.6 s
% Wiki: Interpolate: H =2 m, T = 6.3 s

A_w = 0*1.5/2;          % wave height (m)
w_dir = -1;           % wave direction (1 for 0 deg, -1 for 180 deg)
H = 10;               % mean water level (m)
T_w = 5.7;            % wave period (s)
omega_w = 2*pi/T_w;   % wave circular frequency (rad/s)
k_w = omega_w^2/g;    % wave number (rad/m)
epsilon_w = 0;        % random phase angle (rad)
C_w = omega_w/k_w;    % wave speed (m/s)
lambda_w = C_w*T_w;   % wave length (m)
%% Water Current:
u_t = 0.2*1;          % tidal current component (m/s)
u_w = 0*0.2;         % local wind current component (m/s)

%% floating buoy: two cylinder shape buoys rigidly connected
% the 2 cylinder structure maitain the stability of the buoy
% and prevent roll motion
rho_w = 1000;   % water density (kg/m3)

% Cylinder:
% R_b = 0.1;             % buoy radius (m), cylinder
% A_b = pi*R_b^2;        % (m^2)
% l_b = 0.6;             % buoy length (m)
% m_b = 5;               % (kg)
% rho_b = m_b/(A_b*l_b); %  buoy density (kg/m3)
% J_b = m_b*R_b^2;

% A_im = (m_b/2)/rho_w/l_b; % this is only for 1 one of 2 cylinders
% theta = fzero(@(theta) R_b^2/2*(theta-sin(theta))-A_im,0);
% zeta_im = R_b*(1-cos(theta/2)); % R-H, water surface to buoy tip
% 
% u_s = A_w^2*omega_w*k_w*exp(2*k_w*(-zeta_im/2)); % stocks drift (m/s)
% 
% X0_b = [0,H+R_b-zeta_im];
% V0_b = [omega_w*A_w*cos(k_w*X0_b(1)+epsilon_w), 0];

% Cuboid:
h_b =0.25;             % buoy side (m), cuboid 
A_b = h_b^2;           % (m^2)
l_b = 0.8;             % buoy length (m)
Vol_b = A_b*l_b;       % buoy volume (m)
% for the buoy to be half immersed:
m_b = rho_w*Vol_b/4;   % buoy mass (kg)
J_b = m_b*h_b*l_b/12; 

A_im = (m_b)/rho_w/l_b; % immersed frontal area (m^2) (different from whetted area)
Delta_h = (0.5-A_im/A_b)*h_b; % z_b-zeta, buoy center above water surface 

u_s = A_w^2*omega_w*k_w*exp(2*k_w*(-Delta_h/2)); % stocks drift (m/s)

X0_b = [0,H+Delta_h];
v_x_w = w_dir*omega_w*A_w*sin(k_w*X0_b(1)+epsilon_w);
v_z_w = w_dir*omega_w*A_w*cos(k_w*X0_b(1)+epsilon_w); % ???stocks drift (m/s)
V0_b = [u_t+u_w+u_s+v_x_w, v_z_w];

%Cd_f =3; % drag
%drag_wave = 0.3; % wave lateral drag force (N)

a_11 = 0.05*m_b;        % added mass
a_33 = m_b;
b_11 = 0;               % added dampind
omega_h = omega_w;
b_33 = 2*m_b*omega_h;   

C_S = [5;7]*10^-3;       % skin friction constant
D_s = [2;0.3*3];         % sfkin friction coefficient (viscous)
D_x = b_11 + D_s(1);     % total drag coefficient
D_z = b_33 + D_s(2);

%% Cable:
Lc_0 = 7;                % Cable free length (m)
dL_max = 0.0;            % maximum cable elongation (m)
T_max = 80;              % max cable tension (N)
Kc = T_max/dL_max;       % Cable spring constant (N/m)

alpha_bar_0 = 45*(pi/180);
%zq_bar = X0_b(2) + Lc_0*tan(alpha_bar_0)*sqrt(1-cos(alpha_bar_0)^2)
zq_bar = X0_b(2) + Lc_0*sin(alpha_bar_0);

%% Tension and azimuth control:
kp_a = 10;
kd_a = 10;
ki_a = 2;
kp_T = 4;
kd_T = 0.1;
ki_T = 0.4;

%% quadrotor:

J_u = 0.03;              % moment of inertia (kg.m2)
m_u = 1.8  ;             % quadcopter mass (kg)  1.634
K = 20;                  % motors maximum thrust,each (N)
L = 0.2;                 % qadrotor arm length (m)
Theta0 = atan2(V_bar*D_x*cos(alpha_bar_0),V_bar*D_x*sin(alpha_bar_0)+m_u*g*cos(alpha_bar_0));%10*pi/180;

Tm = 1/15;

%X0_q = X0_b+ [0,3];    % initial CG position (m)
%V0_q = [0,0];          % initial CG velocity (m/s) 

% Limits
Xu_min= [-30 H];        % Window dimensions: min borders
Xu_max= [30 20];        % Window dimensions: max borders
LIMIT_CMD_PITCH = 45*pi/180;
LIMIT_uCMD_HEIGHT = 3.5*K;
LIMIT_CMD_x = 3;

% Aerodynamics:
rho_a = 1.225;             % air density (kg/m^3)
Cd_u = 0.05;               % drag coefficient % 1.05
Ac = [0.0331 0.0331 0.05]; % Cross section area of the quadrotor in x,y and z direction (m^2)

% Gains:
% X:  
%k1 = 3.3; 
%k2 = 0.58; 
kp_t = 10; %k1 = 0.3
kd_t = 6;%k2 = 30
gamma_d_x=0.4;
dc_M_x=0.3;
% Pitch: 
k3 = 4.5;%30;
k4 = 2;% 0.3;  
gamma_d_p=2*0.3/0.03*2*0.2*12/4; %k1=14.6
dc_M_p=6.7;

%motors discrete transfer function:
sys=tf(1,[Tm 1]);
sysd = c2d(sys,T);
Am = sysd.denominator{1};
Bm = sysd.numerator{1};

% Low pass filter:

A_1df = [1/40 1];                   % filter for 1 derivative
A_2df = conv([1/40 1],[1/60 1]);    % filter for 2 derivatives

sys=tf(1,A_1df);
sysd = c2d(sys,T);
Am1_LPF = sysd.denominator{1};
Bm1_LPF = sysd.numerator{1};

sys=tf(1,A_2df);
sysd = c2d(sys,T);
Am2_LPF = sysd.denominator{1};
Bm2_LPF = sysd.numerator{1};

% Z filter
sys=tf(1,conv([1/20 1],[1/20 1]));
sysd = c2d(sys,T);
AmZ_LPF = sysd.denominator{1};
BmZ_LPF = sysd.numerator{1};

% V_bar filter
sys=tf(1,conv([1/0.5 1],[1/0.5 1]));
sysd = c2d(sys,T);
AmV_LPF = sysd.denominator{1};
BmV_LPF = sysd.numerator{1};

% alpha_dot filter
sys=tf(1,conv([1/5 1],[1/5 1]));
sysd = c2d(sys,T);
AmAlpha_LPF = sysd.denominator{1};
BmAlpha_LPF = sysd.numerator{1};

%
epsilon_2 = 0.05*Vol_b;
A_whetted =  l_b*h_b + 2*l_b*(h_b/2) % average
D_from_C = A_whetted*V_bar*0.5*rho_w*C_S(1)
v_x_w_max = omega_w*A_w;

V_min = epsilon_1*cosd(alpha_bar)/D_from_C + u_t + v_x_w_max
V_max = (m_b + m_u - epsilon_2*rho_w)*g*tan(Theta0)...
            /D_x+u_t+u_w+u_s+abs(v_x_w)
