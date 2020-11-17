function psi_alpha = finite_time_feedback_law(alpha)

phi = @(x1, x2) x1 + (1/2 - alpha)*sign(x2)*abs(x2)^(2 - alpha);
psi_alpha = @(x1, x2) -sign(x2)*abs(x2)^alpha ...
                    - sign(phi(x1, x2))*abs(phi(x1, x2))^(-alpha/2);
    