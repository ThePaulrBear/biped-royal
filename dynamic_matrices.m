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
s12 = sin(theta1 - theta2);
s13 = sin(theta1 - theta3);
c12 = cos(theta1 - theta2);
c13 = cos(theta1 - theta3);
L = torso_length; % this is a lowercase L in Grizzle.
r = leg_length;
D = [(5/4*m + MH + MT)*r^2, -1/2*m*r^2*c12, MT*r*L*c13;
            -1/2*m*r^2*c12,      1/4*m*r^2,          0;
                MT*r*L*c13,              0,     MT*L^2];

C = [                   0, -1/2*m*r^2*s12*omega2, MT*r*L*s13*omega3;
     1/2*m*r^2*s12*omega1,                     0,                 0;
       -MT*r*L*s13*omega1,                     0,                 0];

G = [-1/2*g*(2*MH + 3*m+ 2*MT)*r*sin(theta1);
                      1/2*g*m*r*sin(theta2);
                        -g*MT*L*sin(theta3)];