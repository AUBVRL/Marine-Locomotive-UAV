% matrix calculation

syms L_f L_d
% e_L = L_d - L_f
G0=[-cos(L_d - L_f) -sin(L_d - L_f)
     sin(L_d - L_f) -cos(L_d - L_f)];
 
 
G1=[ cos(L_d) -sin(L_d)
     sin(L_d)  cos(L_d)];
 
 
G = simplify(inv(G0*G1))
