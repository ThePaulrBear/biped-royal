% function u = control_law(x)
% 
% % function [C, D, g_sig] = dynamic_matrices(x)
% 
%     global g_sig torso_mass_sig hip_mass_sig leg_mass_sig leg_length_sig torso_length_sig
%     theta1 = x(1);
%     theta2 = x(2);
%     theta3 = x(3);
%     omega1 = x(4);
%     omega2 = x(5);
%     omega3 = x(6);
% 
%     MT = 1;
%     MH = 1;
%     m = 1;
%     
%     s12 = sin(theta1 - theta2);
%     s13 = sin(theta1 - theta3);
%     c12 = cos(theta1 - theta2);
%     c13 = cos(theta1 - theta3);
%     g_sig = 1;
%     L = 1; % this is a lowercase L in Grizzle.
%     r = 1;
%     D = [(5/4*m + MH + MT)*r^2, -1/2*m*r^2*c12, MT*r*L*c13;
%                 -1/2*m*r^2*c12,      1/4*m*r^2,          0;
%                     MT*r*L*c13,              0,     MT*L^2];
% 
%     C = [                   0, -1/2*m*r^2*s12*omega2, MT*r*L*s13*omega3;
%          1/2*m*r^2*s12*omega1,                     0,                 0;
%            -MT*r*L*s13*omega1,                     0,                 0];
% 
%     G = [-1/2*g_sig*(2*MH + 3*m+ 2*MT)*r*sin(theta1);
%                           1/2*g_sig*m*r*sin(theta2);
%                             -g_sig*MT*L*sin(theta3)];
% % end
% 
% global B_sig
% omega = x(4:6);
% omega2 = omega(2);
% omega3 = omega(3);
% 
% A = [0, 0, 1; 
%      0, 1, 1];
% % [C, D, g_sig] = dynamic_matrices(x);
% 
% D_inv = inv(D);
% Lf = [omega3; omega2 + omega3]; 
% LfLf = A * D_inv *(-C*omega - G);
% LgLf = A * D_inv * B_sig;
% 
% theta1 = x(1);
% theta2 = x(2);
% theta3 = x(3);
% % [Lf, LfLf, LgLf] = Lie_derivatives(x);
% 
% y = [theta3 + pi / 2; theta1 + theta2];
% ydot = Lf;
% 
% alpha = 0.5;
% phi_alpha = y + (1/2 - alpha)*sign(ydot).*abs(ydot).^(2 - alpha);
% Psi = -sign(ydot).*abs(ydot).^alpha - sign(phi_alpha).*abs(phi_alpha).^(alpha/(2-alpha));
% 
% % Psi(isnan(Psi)) = 0;
% 
% u = -inv(LgLf)\(Psi - LfLf);
% 
% end
% 
% function y = phi(x)
%     y = x(1) + (1/2 - alpha)*sign(x(2)).*abs(x(2)).^(2 - alpha);
% end
% 
% function y = psi_f(x)
%      y = -sign(x(2)).*abs(x(2)).^alpha - sign(phi(x)).*abs(phi(x)).^(alpha/(2-alpha));
% end