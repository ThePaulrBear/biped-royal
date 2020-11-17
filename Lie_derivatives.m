function [Lf, LfLf, LgLf] = Lie_derivatives(x)

global B
omega = x(4:6);
omega2 = x(2);
omega3 = x(3);

A = [0, 0, 1; 
     0, 1, 1];
[C, D, G] = dynamic_matrices(x);
 
D_inv = inv(D);
Lf = [omega3; omega2 + omega3]; 
LfLf = A * D_inv *(-C*omega - G);
LgLf = A * D_inv * B;