

%% cylinder
omega = 2*pi/15;
rho_w = 1000;
R_b = 0.1;
l_b = 0.6;
V_b = 2*pi*R_b^2*l_b;
m_b = 5;
rho_b = m_b/V_b;
A_im = m_b/rho_w/l_b;
% b_33 = 2*rho_w*A_im*l_b*omega;
b_33 = 2*m_b*omega

theta = fzero(@(theta) R_b^2/2*(theta-sin(theta))-A_im,0);
zeta_im = R_b*(1-cos(theta/2)) % R-H, water surface to buoy tip

%% Skin friction
h_b =0.15;             % buoy side (m), cuboid 
l_b = 0.5;             % buoy length (m)
nu_w = 1.787*10^-6; % kinematic viscosity of water
V = [0.1:0.5:5];
for i=1:length(V)
    V_b = V(i); % buoy average velocity (m/s)
    Re = V_b*l_b/nu_w;
    C_S(i,1) = 0.075/(log10(Re)-2)^2 % skin friction coeff
    % F_D(i,1) = C_F*(0.5*rho_w*V_b^2)*(2*R_b*theta*l_b); % cylinder
    F_D(i,1) = C_S(i,1)*(0.5*rho_w*V_b^2)*(l_b*h_b + l_b*h_b/2 + l_b*h_b/2); % Cuboid
end

figure(1)
plot(V,F_D)
figure(2)
plot(V,C_S)


%% Steady state and tension:
V_bar_tilde = (V_bar-u_t-u_s)
T_bar = D_x*V_bar_tilde/cos(alpha_bar)
theta_bar = atan(D_x*V_bar_tilde*cos(alpha_bar)/(m_u*g*cos(alpha_bar)+D_x*V_bar_tilde*sin(alpha_bar)))
theta_bar = atan(D_x*V_bar_tilde*sin(alpha_bar)/(m_u*g*sin(alpha_bar)-D_x*V_bar_tilde*cos(alpha_bar)+T_bar))

theta_bar_d = theta_bar*180/pi
u_1 = m_u*g*cos(alpha_bar)/cos(alpha_bar+theta_bar)
u_1 = D_x*V_bar_tilde/sin(theta_bar)
u_1 = (T_bar+m_u*g*sin(alpha_bar))/sin(alpha_bar+theta_bar)
% function F = root2d(x)
% 
% F(1) = D_11_I*x_dot_bar-x(1)*sin(x(2));
% F(2) = D_21_I*x_dot_bar+(m_q+m_b)*g-x(1)*cos(x(3))-rho_w*x(3)*g;
% F(3) = m_q*g*Lc_0*cos(alpha_bar)-x(1)*Lc_0*cos(alpha_bar+x(2))
% 
% fun = @root2d;
% x0 = [0,0];
% x = fsolve(fun,x0,options)

