close all
% clc 
startup;
global alpha theta3_ref
global x0
alpha = 0.5;
theta3_ref = 0.0;

simulation = "continuous";

switch simulation
    case "test"
        disp("hi!")
    case "1D"
        x0 = ones(2, 1);
        [t, y] = ode45(@f_linear, [0, 5], x0);
        plot(t, y);
    case "continuous"
        [t, x] = ode23(@f_cl, [0, 5], x0);

        u = [0*t, 0*t];
        for i = 1:length(t)
            u(i, :) = control_law(x(i, :)')';
        end
        
        fprintf("\n")
        fprintf("\\theta_{10} = %.1f, \\theta_{20} = %.1f, \\theta_{30} = %.1f\n", ...
                    x(1,1), x(1,2), x(1,3));
        fprintf("r = %.1f, \\ell = %.1f\n", leg_length, torso_length);
        
        plot(t, u)
        title("Control u(t)")
        legend("u_1", "u_2")
        xlabel("t")
        
        figure()
        y = [x(:, 3), x(:, 1) + x(:, 2)];
        plot(t, y);
        title("Output y(t)")
        legend("\theta_3", "\theta_1 + \theta_2")
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
    u = control_law(x);
    omega = x(4:6);
    xdot = [omega; 
            D \ (-C*omega - G + B*u)]; 
end

function u = control_law(x)
 
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

    
    % Center of gravity constraint from page 7 of Grizzle
    assert(0 < L * MT)
    assert(L*MT < r*(m + MT + MH))

    A = [0, 0, 1; 
         1, 1, 0];
    global B
    [C, D, G] = dynamic_matrices(x);

    Lf = [omega3; omega1 + omega2]; 
    LfLf = A * (D \ (-C*omega - G));
    LgLf = A * (D \ B);

    R11 = m*r^3/4*(5/4*m*r + MH*r + MT*r - m*r*c12^2 + MT*L*c13);
    R12 = m*r^3/4*(5/4*m*r + MH*r + MT*r - m*r*c12^2 + 2*MT*L*c12*c13);
    R21 = -m*MT*L*r^2/4*(1 + 2*c12)*(r*c13 + L);
    R22 = -MT*L*r^2 / 4 *(5*m*L + 4*MH*L + 4*MT*L + m*r*c13 + 2*m*r*c12*c13 - 4*MT*L*(c13)^2 + 2* m*L*c12);
    detD = m* MT*r^4*L^2/4*(5/4*m + MH + MT - m*c12^2 - MT*c13^2);
%     assert(abs(detD - det(D)) < 1e-6)
    LgLf_alt = 1/detD * [R11, R12; R21, R22];
%     assert(norm(LgLf - LgLf_alt) < 1e-6)

    y = [theta3 - theta3_ref; theta1 + theta2];
    ydot = Lf;

    Psi = [psi_f([y(1); ydot(1)]); 
           psi_f([y(2); ydot(2)])];
    u = LgLf \ (Psi - LfLf);

end

function y = phi(x)
    global alpha
    y = x(1) + (1/2 - alpha)*sign(x(2)).*abs(x(2)).^(2 - alpha);
end

function y = psi_f(x)
    global alpha
     y = -sign(x(2)).*abs(x(2)).^alpha - sign(phi(x)).*abs(phi(x)).^(alpha/(2-alpha));
end

function [C, D, G] = dynamic_matrices(x)
    global g torso_mass hip_mass leg_mass leg_length torso_length
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    s12 = sin(theta1 + theta2);
    s13 = sin(theta1 + theta3);
    c12 = cos(theta1 - theta2);
    c13 = cos(theta1 - theta3);
    
    MT = torso_mass;
    MH = hip_mass;
    m = leg_mass;
    L = torso_length; % this is a lowercase L in Grizzle.
    r = leg_length;

    D = [(5*m/4 + MH + MT)*r^2, -1/2*m*r^2*c12, MT*r*L*c13;
                -1/2*m*r^2*c12,      1/4*m*r^2,          0;
                    MT*r*L*c13,              0,     MT*L^2];

    C = [                   0, -1/2*m*r^2*s12*omega2, MT*r*L*s13*omega3;
         1/2*m*r^2*s12*omega1,                     0,                 0;
           -MT*r*L*s13*omega1,                     0,                 0];

    G = g*[-1/2*(2*MH + 3*m + 2*MT)*r*sin(theta1);
                              1/2*m*r*sin(theta2);
                                -MT*L*sin(theta3)];
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

