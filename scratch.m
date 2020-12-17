close all
% clc 
startup;
global alpha theta3_ref
global x0
alpha = 0.5;
theta3_ref = 7*pi/8;

simulation = "continuous";

switch simulation
    case "test"
        disp("hi!")
    case "1D"
        x0 = ones(2, 1);
        [t, y] = ode45(@f_linear, [0, 5], x0);
        plot(t, y);
    case "continuous"
        [t, x] = ode23(@f_cl, [0, 10], x0);

        u = 0*t*control_law(x(1, :)')';
        for i = 1:length(t)
            u(i, :) = control_law(x(i, :)')';
        end
        
        fprintf("\n")
        fprintf("\\theta_{10} = %.1f, \\theta_{20} = %.1f, \\theta_{30} = %.1f\n", ...
                    x(1,1), x(1,2), x(1,3));
        fprintf("r = %.1f, \\ell = %.1f\n", leg_length, torso_length);
        
        figure
        plot(t, u)
        title("Control u(t)")
        legend("\tau_1", "\tau_2", "\tau_3", "\tau_4", "\tau_5")
        xlabel("t")
        
        figure
        subplot(2, 1, 1)
        plot(t, x(:,1:5))
        title("Angles \theta(t)")
        legend("\theta_1", "\theta_2", "\theta_3", "\theta_4", "\theta_5")
        xlabel("t")
        subplot(2, 1, 2)
        plot(t, x(:,6:10))
        title("Anglar Acceleration \omega(t)")
        legend("\omega_1", "\omega_2", "\omega_3", "\omega_4", "\omega_5")
        xlabel("t")
        
        figure()
        y = [x(:, 1) - pi/2, x(:, 2), x(:, 3) - theta3_ref, x(:, 5)];
        plot(t, y);
        title("Output y(t)")
        legend("\theta_1 - \pi/2", "\theta_2", "\theta_3 - \theta_3^d", "\theta_4", "\theta_5")
        xlabel("t")
        ylabel("y")
        
%         figure()
%         subplot(2, 1, 1)
%         plot(t, x(:, 1:3))
%         title("State x(t)")
%         legend("\theta_1", "\theta_2", "\theta_3")
%         xlabel("t")
%         
%         subplot(2, 1, 2)
%         plot(t, x(:, 4:6))
%         legend("\omega_1", "\omega_2", "\omega_3")
%         xlabel("t")
    case "hybrid"
        % initial conditions;
        % simulation horizon
        TSPAN=[0 10];
        JSPAN = [0 20];
        rule = 1;
        options = odeset('RelTol',1e-6,'MaxStep',.1);
        % simulate
        [t,j,x] = HyEQsolver(@flow_map,@jump_map,@C,@D,x0,TSPAN,JSPAN,rule,options);
        % plot solution
        figure(1) % position
        clf
        subplot(2,1,1),plotflows(t,j,x(:,1))
        grid on
        ylabel('\theta_1')
        subplot(2,1,2),plotflows(t,j,x(:,2))
        grid on
        ylabel('\theta_2')
        figure(2) % velocity
        % plot hybrid arc
        plotHybridArc(t,j,x)
        xlabel('j')
        ylabel('t')
        zlabel('x1')
end

function xdot = f_linear(~, x)
    xdot = [x(2); psi_f(x)]; 
end

function xdot = f_cl(~, x)
    global B
    [C, D, G] = dynamic_matrices(x);
    [u, v] = control_law(x);
    omega = x(6:10);
    xdot = [omega; 
            D \ (-C*omega - G + B*u)]; 
%     assert(norm(xdot(6) - v(1)) < 0.01);


    % Check our calculations.
    omega_dot = xdot(6:10);
    
    global g torso_mass leg_mass leg_length torso_length
    m1 = leg_mass; %m1_sig
    m2 = leg_mass; %m2_sig
    m3 = leg_mass; %m3_sig
    m4 = leg_mass; %m4_sig
    m5 = torso_mass; %m5_sig
    
    L1 = leg_length; %L1_sig
    L2 = leg_length; %L2_sig
    L3 = leg_length; %L3_sig
    L4 = leg_length; %L4_sig
    L5 = torso_length; %L5_sig   
    
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    theta4 = x(4);
    theta5 = x(5);
    
    % alphas
    alpha10 = theta1;
    
    alpha21 = theta2;
    alpha31 = alpha21+theta3;
    alpha41 = alpha31+theta4;
    alpha51 = alpha21+theta5;
    
    alpha32 = theta3;
    alpha42 = alpha32+theta4;
    alpha52 = theta5;
    
    alpha43 = theta4;
    
    alpha34 = -alpha43;
    
    alpha25 = -alpha52;
    
    R21 = rotation_2D(alpha21);
    R31 = rotation_2D(alpha31);
    R41 = rotation_2D(alpha41);
    R51 = rotation_2D(alpha51);
    
    R32 = rotation_2D(alpha32);
    R42 = rotation_2D(alpha42);
    R52 = rotation_2D(alpha52);
    
    R12 = R21';
    R23 = R32';
    R34 = rotation_2D(alpha34);
    R25 = rotation_2D(alpha25);
    
    R43 = rotation_2D(alpha43);
    
    O = zeros(2);
    
    A = [O, O, O, O, O;
         R21, O, O, O, O;
         R31, R32, O, O, O;
         R41, R42, R43, O, O;
         R51, R52, O, O, O];
    P = [O, R12, O, O, O;
         O, O, R23, O, R25;
         O, O, O, R34, O;
         O, O, O, O, O;
         O, O, O, O, O];
    L = diag([L1, L1, L2, L2, L3, L3, L4, L4, 2*L5, 2*L5]);
    M = diag([m1, m1, m2, m2, m3, m3, m4, m4, m5, m5]);
    E2 = zeros(5, 10);
    for i = 1:5
        row = i;
        col = 2*i;
       E2(row, col) = 1;
    end
    w = NaN*zeros(10, 1);
    for i = 1:5
        w(2*i - 1) = -omega(i)^2;
        w(2*i) = omega_dot(i);
    end
    v_c_dot = (A + 1/2*eye(10)) * L * w;
    f = inv(eye(10) - P) * M * v_c_dot;
    torque_eq = 1/2 * E2 * M * L * v_c_dot + E2 * L * P * f + G - B * u;
%     eig(D)
%     eig(C)
    
    assert(norm(u) < 1e6, "control is very large!")
    assert(norm(torque_eq) < 0.01, "Torque equation is incorrect");

end

function [u, v] = control_law(x)
 
    global theta3_ref B

    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    theta4 = x(4);
    theta5 = x(5);
    omega = x(6:10);

    y = [theta1 - pi/2;
         theta2;
         theta3 - theta3_ref;
         theta4;
         theta5];
    H = [1 0 0 0 0; 
         0 1 0 0 0; 
         0 0 1 0 0; 
         0 0 0 1 0; 
         0 0 0 0 1];

    %hip_mass_sig is removed from the parameters of dynamic_matrices
    [C, D, G] = dynamic_matrices(x);

    Lf = H*omega;
    LfLf = H * (D \ (-C*omega - G));
    LgLf = H * (D \ B);

    ydot = Lf;

    global alpha
    v = [psi_f([y(1); ydot(1)], alpha); 
         psi_f([y(2); ydot(2)], alpha); 
         psi_f([y(3); ydot(3)], alpha); 
         psi_f([y(4); ydot(4)], alpha); 
         psi_f([y(5); ydot(5)], alpha)
         ];
    u = LgLf \ (v - LfLf);
end

function y = phi(x)
     global alpha
    y = x(1) + (1/2 - alpha)*sign(x(2)).*abs(x(2)).^(2 - alpha);
end

function y = psi_f(x, alpha)
     y = -sign(x(2)).*abs(x(2)).^alpha - sign(phi(x)).*abs(phi(x)).^(alpha/(2-alpha));
end

function [C, D, G] = dynamic_matrices(x)
                    %m1_sig, m2_sig, m3_sig, m4_sig, m5_sig, L1_sig, L2_sig, L3_sig, L4_sig, L5_sig
    global g torso_mass leg_mass leg_length torso_length

    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    theta4 = x(4);
    theta5 = x(5);
    omega = x(6:10);
    omega1 = x(6);
    omega2 = x(7);
    omega3 = x(8);
    omega4 = x(9);
    omega5 = x(10);
    
    % alphas
    alpha10 = theta1;
    alpha20 = alpha10+theta2;
    alpha30 = alpha20+theta3;
    alpha40 = alpha30+theta4;
    alpha50 = alpha20+theta5;
    
    alpha21 = theta2;
    alpha31 = alpha21+theta3;
    alpha41 = alpha31+theta4;
    alpha51 = alpha21+theta5;
    
    alpha12 = -alpha21;
    alpha32 = theta3;
    alpha42 = alpha32+theta4;
    alpha52 = theta5;
    
    alpha13 = -alpha31;
    alpha23 = -alpha32;
    alpha43 = theta4;
    
    alpha14 = -alpha41;
    alpha24 = -alpha42;
    alpha34 = -alpha43;
    
    alpha15 = -alpha51;
    alpha25 = -alpha52;
    
    %R_22 = cos()
    R10_22 = cos(alpha10);
    R20_22 = cos(alpha20); 
    R30_22 = cos(alpha30); 
    R40_22 = cos(alpha40);
    R50_22 = cos(alpha50);
    
    R21_22 = cos(alpha21);
    R31_22 = cos(alpha31);
    R41_22 = cos(alpha41);
    R51_22 = cos(alpha51);
    
    R12_22 = cos(alpha12);
    R32_22 = cos(alpha32);
    R42_22 = cos(alpha42);
    R52_22 = cos(alpha52);
    
    R13_22 = cos(alpha13);
    R23_22 = cos(alpha23);
    R43_22 = cos(alpha43);
    
    R14_22 = cos(alpha14);
    R24_22 = cos(alpha24);
    R34_22 = cos(alpha34);
    
    R15_22 = cos(alpha15);
    R25_22 = cos(alpha25);
    
    %R_21 = sin()
    R21_21 = sin(alpha21);
    R31_21 = sin(alpha31);
    R41_21 = sin(alpha41);
    R51_21 = sin(alpha51);
    
    R12_21 = sin(alpha12);
    R32_21 = sin(alpha32);
    R42_21 = sin(alpha42);
    R52_21 = sin(alpha52);
    
    R13_21 = sin(alpha13);
    R23_21 = sin(alpha23);
    R43_21 = sin(alpha43);
    
    R14_21 = sin(alpha14);
    R24_21 = sin(alpha24);
    R34_21 = sin(alpha34);
    
    R15_21 = sin(alpha15);
    R25_21 = sin(alpha25);
    
    R10 = rotation_2D(alpha10);
    R20 = rotation_2D(alpha20);
    R30 = rotation_2D(alpha30);
    R40 = rotation_2D(alpha40);
    R50 = rotation_2D(alpha50);
    
    R21 = rotation_2D(alpha21);
    R31 = rotation_2D(alpha31);
    R41 = rotation_2D(alpha41);
    R51 = rotation_2D(alpha51);
    
    R32 = rotation_2D(alpha32);
    R42 = rotation_2D(alpha42);
    R52 = rotation_2D(alpha52);
    
    R12 = R21';
    R23 = R32';
    R34 = rotation_2D(alpha34);
    R25 = rotation_2D(alpha25);
    
    R43 = rotation_2D(alpha43);
    
    %I defined the masses and lengths in terms of leg_mass and
    %leg_length because they're already defined
    m1 = leg_mass; %m1_sig
    m2 = leg_mass; %m2_sig
    m3 = leg_mass; %m3_sig
    m4 = leg_mass; %m4_sig
    m5 = torso_mass; %m5_sig
    
    L1 = leg_length; %L1_sig
    L2 = leg_length; %L2_sig
    L3 = leg_length; %L3_sig
    L4 = leg_length; %L4_sig
    L5 = torso_length; %L5_sig   
    
    O = zeros(2);
    
    A = [O, O, O, O, O;
         R21, O, O, O, O;
         R31, R32, O, O, O;
         R41, R42, R43, O, O;
         R51, R52, O, O, O];
    P = [O, R12, O, O, O;
         O, O, R23, O, R25;
         O, O, O, R34, O;
         O, O, O, O, O;
         O, O, O, O, O];
    M = diag([m1, m2, m3, m4, m5]);
    M_10x10 = diag([m1, m1, m2, m2, m3, m3, m4, m4, m5, m5]);
    L = diag([L1, L2, L3, L4, 2*L5]);
    L_10x10 = diag([L(1,1), L(1,1), L(2,2), L(2,2), L(3, 3), L(3, 3), L(4, 4), L(4, 4), L(5, 5), L(5, 5)]);
    
    E2 = zeros(5, 10);
    for i = 1:5
        row = i;
        col = 2*i;
       E2(row, col) = 1;
    end
    
    Pi = [(m1/4)+m2+m3+m4+m5, (m2+m3+m4+m5)*R12_22, (m3+m4)*R13_22, m4*R14_22, m5*R15_22;
        ((m2/2)+m3+m4+m5)*R21_22, m2/4+m3+m4+m5, (m3+m4)*R23_22, m4*R24_22, m5*R25_22;
        ((m3/2)+m4)*R31_22, ((m3/2)+m4)*R32_22, (m3/4)+m4, m4*R34_22, 0;
        m4/2*R41_22, m4/2*R42_22, m4/2*R43_22, m4/4, 0;
        m5/2*R51_22, m5/2*R52_22, 0, 0, m5/2];
    
    Lambda = [0, (m2+m3+m4+m5)*R12_21, (m3+m4)*R13_21, m4*R14_21, m5*R15_21;
        ((m2/2)+m3+m4+m5)*R21_21, 0, (m3+m4)*R23_21, m4*R24_21, m5*R25_21;
        ((m3/2)+m4)*R31_21, ((m3/2)+m4)*R32_21, 0, m4*R34_21, 0;
        m4/2*R41_21, m4/2*R42_21, m4/2*R43_21, 0, 0;
        m5/2*R51_21, m5/2*R52_21, 0, 0, 0];
    
    G = g*[L1*(3/2*m1+m2+m3+m4+m5)*R10_22;
           L2*(3/2*m2+m3+m4+m5)*R20_22;
           L3*(3/2*m3+m4)*R30_22;
           L4*3/2*m4*R40_22;
           2*L5*m4*R50_22];       
    Gamma = [R10; R20; R30; R40; R50]*[0; g];
    G_other = (1/2 * M * L * E2 + L * E2 * P * inv(eye(10) - P) * M_10x10) * Gamma;
    G = G_other;
%     assert(norm(G - G_other) < 0.1, "G's don't match!")
       
    D = L*Pi*L;    
    C = -L*Lambda*L*diag([omega1, omega2, omega3, omega4, omega5]);
    
    some_5x10 = 1/2*M*L*E2*(A+1/2*eye(10))*L_10x10 + L*E2*P*inv(eye(10) - P)*M_10x10*(A + 1/2*eye(10))*L_10x10;
    D_other = some_5x10(:, 2:2:10);
    C_other = -some_5x10(:, 1:2:10)*diag(omega);
%     assert(norm(C_other * omega + D_other * omega_dot - some5x10*w) < 1e-6)
%     D - D_other;
    C = C_other;
    D = D_other;
    
    assert(all(eig(D) > 0), "D is not positive definite!")
end

 function xdot = flow_map(x)
    xdot = f_cl(NaN, x);
 end
 
function xplus = jump_map(x)
    global torso_mass hip_mass leg_mass leg_length torso_length theta3_ref

    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega = x(4:6);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);
    
    c12 = cos(theta1 - theta2);
    c13 = cos(theta1 - theta3);

    MT = torso_mass;
    MH = hip_mass;
    m = leg_mass;
    L = torso_length; % this is a lowercase L in Grizzle.
    r = leg_length;
    
    theta1 = x(2);
    theta2 = x(1);
    theta3 = x(3);
    
    De = NaN*zeros(5);
    De(1,1) =((5)/(4)*m+MH+MT)*r^(2);
    De(1,2) = -(1)/(2)*m*r^(2)*c12;
    De(2,1) = De(1,2);
    De(1,3) = MT*r*L*c13;
    De(3,1) = De(1,3);
    De(1,4) =((3)/(2)*m+MH+MT)*r*cos(theta1);
    De(4,1) = De(1,4);
    De(1,5) =-((3)/(2)*m+MH+MT)*r*sin(theta1);
    De(5,1) = De(1,5);
    De(2,2) =(1)/(4)*m*r^(2);
    De(2,3) =0;
    De(3,2) = De(2,3);
    De(2,4) =-(1)/(2)*m*r*cos(theta2);
    De(4,2) = De(2,4);
    De(2,5) = (1)/(2)*m*r*sin(theta2);
    De(5,2) = De(2,5);
    De(3,3) = MT*L^(2);
    De(3,4) = MT*L*cos(theta3);
    De(4,3) = De(3,4);
    De(3,5) = -MT*L*sin(theta3);
    De(5,3) = De(3,5);
    De(4,4) = 2*m+MH+MT;
    De(4,5) = 0;
    De(5,4) = De(4,5);
    De(5,5) = 2*m+MH+MT;
    assert(norm(De - De') < 1e-10)
    
    E = [r*cos(theta1), -r*cos(theta2), 0, 1, 0;
        -r*sin(theta1),  r*sin(theta2), 0, 0, 1];
    xplus = zeros(size(x));
    xplus(1) = x(2);
    xplus(2) = x(1);
    xplus(3) = x(3);
    xplus(4) = x(5);
    xplus(5) = x(4);
    xplus(6) = x(6);
    
end
 
function is_flow = C(x)
    theta1 = x(1);
    theta2 = x(2);
    
    is_flow = (theta2 > theta1); 
end

function is_jump = D(x)
    theta1 = x(1);
    theta2 = x(2);
    
    is_jump = (theta2 <= theta1); 
end

function R = rotation_2D(alpha)
    R = [cos(alpha), -sin(alpha);
         sin(alpha), cos(alpha)];
end
