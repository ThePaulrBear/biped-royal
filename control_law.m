function u = control_law(x)

[Lf, LfLf, LgLf] = Lie_derivatives(x);


u = inv(LgLf)*(Psi - LfLf)