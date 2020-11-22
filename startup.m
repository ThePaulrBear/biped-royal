global g torso_mass hip_mass leg_mass leg_length torso_length B

g = -9.8; % Check the sign for this!
torso_mass = 0.3; % "M_T"
hip_mass = 2.0; % "M_H"
leg_mass = 1.0; % "m" (concentrated at center of mass))

leg_length = 1.0;
leg_radius = 0.01;
torso_length = 0.5;
% torso_radius = 0.2; % not implemented?

B = [-1,  0;
      0, -1; 
      1, 1];

alpha = 0.5;