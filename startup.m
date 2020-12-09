global g torso_mass hip_mass leg_mass leg_length torso_length B 
global theta_A_initial theta_B_initial torso_initial x0 theta3_ref
global alpha z_start
global static_friction kinetic_friction 

g = 9.80665; % Check the sign for this!
torso_mass = 1.0; % "M_T"
hip_mass = 1.0; % "M_H"
leg_mass = 1.0; % "m" (concentrated at center of mass))

leg_length = 1;
leg_radius = 0.01;
torso_length = 0.6;
% torso_radius = 0.2; % not implemented?

theta_A_initial = 0.1; 
theta_B_initial = -0.6;
theta_Torso_initial = 0.5;
theta3_ref = 0.0;
x0 = [theta_A_initial, theta_B_initial, theta_Torso_initial, 0, 0, 0];

B = [-1,  0;
      0, -1; 
      1, 1];

alpha = 0.9;
z_start = leg_length*max(cos(theta_A_initial), cos(theta_B_initial))+leg_radius;

static_friction = 2.0;
kinetic_friction = 1.0;