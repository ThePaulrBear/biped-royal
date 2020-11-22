close all
clc 
global alpha
alpha = 0.1;

% x0 = ones(2, 1);
% [t, y] = ode45(@f_linear, [0, 5], x0);
% plot(t, y);
% function xdot = f_linear(~, x)
%     xdot = [x(2); psi_f(x)]; 
% end
% 
x0 = 1*rand(6, 1);
[t, x] = ode23(@f_cl, [0, 10], x0);

y = [x(:, 3), x(:, 1) + x(:, 2)];
plot(t, y);
title("output")
figure()
plot(t, x)
title("state")
legend("\theta_1", "\theta_2", "\theta_3", "\omega_1", "\omega_2", "\omega_3")

function xdot = f_cl(~, x)
    global B
    [C, D, G] = dynamic_matrices(x);
    u = control_law(x);
    omega = x(4:6);
    xdot = [omega; 
            inv(D)*(-C*omega - G + B*u)]; 
end

function [C, D, G] = dynamic_matrices(x)
    global g torso_mass hip_mass leg_mass leg_length torso_length
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    MT = torso_mass;
    MH = hip_mass;
    m = leg_mass;
    L = torso_length; % this is a lowercase L in Grizzle.
    r = leg_length;
    
    s12 = sin(theta1 - theta2);
    s13 = sin(theta1 - theta3);
    c12 = cos(theta1 - theta2);
    c13 = cos(theta1 - theta3);

    D = [(5/4*m + MH + MT)*r^2, -1/2*m*r^2*c12, MT*r*L*c13;
                -1/2*m*r^2*c12,      1/4*m*r^2,          0;
                    MT*r*L*c13,              0,     MT*L^2];

    C = [                   0, -1/2*m*r^2*s12*omega2, MT*r*L*s13*omega3;
         1/2*m*r^2*s12*omega1,                     0,                 0;
           -MT*r*L*s13*omega1,                     0,                 0];

    G = g*[-1/2*(2*MH + 3*m+ 2*MT)*r*sin(theta1);
                             1/2*m*r*sin(theta2);
                               -MT*L*sin(theta3)];
end

function u = control_law(x)
 
global g torso_mass hip_mass leg_mass leg_length torso_length
    global B
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    MT = torso_mass;
    MH = hip_mass;
    m = leg_mass;
    L = torso_length; % this is a lowercase L in Grizzle.
    r = leg_length;
    
    s12 = sin(theta1 - theta2);
    s13 = sin(theta1 - theta3);
    c12 = cos(theta1 - theta2);
    c13 = cos(theta1 - theta3);

    [C, D, G] = dynamic_matrices(x);
    omega = x(4:6);

    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    A = [0, 0, 1; 
         1, 1, 0];

    D_inv = inv(D);
    Lf = [omega3; omega1 + omega2]; 
    LfLf = A * (D_inv *(-C*omega - G));
    LgLf = A * (D_inv * B);

    R11 = m*r^3/4*(5/4*m*r + MH*r + MT*r - m*r*c12^2 + MT*L*c13);
    R12 = m*r^3/4*(5/4*m*r + MH*r + MT*r - m*r*c12^2 + 2*MT*L*c12*c13);
    R21 = -m*MT*L*r^2/4*(1 + 2*c12)*(r*c13 + L);
    R22 = -MT*L*r^2 / 4 *(5*m*L + 4*MH*L + 4*MT*L + m*r*c13 + 2*m*r*c12*c13 - 4*MT*L*(c13)^2 + 2* m*L*c12);
    detD = m* MT*r^4*L^2/4*(5/4*m + MH + MT - m*c12^2 - MT*c13^2);
    assert(abs(detD - det(D)) < 1e-6)
    LgLf_alt = 1/detD * [R11, R12; R21, R22];
    assert(norm(LgLf - LgLf_alt) < 1e-6)
    y = [theta3; theta1 + theta2];
    ydot = Lf;

    Psi = [psi_f([y(1); ydot(1)]); psi_f([y(2); ydot(2)])];
    u = inv(LgLf)*(Psi - LfLf);

end

function y = phi(x)
    global alpha
    y = x(1) + (1/2 - alpha)*sign(x(2)).*abs(x(2)).^(2 - alpha);
end

function y = psi_f(x)
    global alpha
     y = -sign(x(2)).*abs(x(2)).^alpha - sign(phi(x)).*abs(phi(x)).^(alpha/(2-alpha));
end