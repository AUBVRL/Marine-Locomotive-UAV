%%

%syms phi(t) theta(t) psi(t) Ixx Iyy Izz
syms phi theta psi phi_d theta_d psi_d Jxx Jyy Jzz
syms p q r T1 T2 T3 Cdx Cdy Cdz x_d y_d z_d v_x v_y v_z
syms dx dy dz ft x y z dp
syms mx mz g x_dot z_dot m11 m22 m12 
J = Jyy;
Theta_d = theta_d;
Omega = q;
T_r = T2;
% Kt = diag([Cdx Cdy Cdz]);
P_d = [x_d; z_d];
V = [v_x; v_z];
Delta = [dx ; dz];
f_r = [0; ft];
X = [x;z;theta];
X_dot = [x_dot;z_dot;q];
M = diag([mx,mz,Jyy]);
G = [0;-g;0];
D = diag([dx;dz;dp]);
R_e2b = [ cos(theta) -sin(theta) 0
          sin(theta)  cos(theta) 0
          0             0        1];
        
R_b2e = [  cos(theta)  sin(theta) 0
          -sin(theta)  cos(theta) 0
          0             0        1]; % J(eta) = J^(-T)
        
D_eta = simplify(simplify(transpose(inv(R_b2e))* D )* simplify(inv(R_b2e)))

M_eta = simplify(R_b2e*M*inv(R_b2e))
M_eta_s = [m11 m12 0; m12 m22 0; 0 0 Jyy];
simplify(simplify(X_dot'*M_eta_s*X_dot))
simplify(simplify(X_dot'*D*X_dot))
