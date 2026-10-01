function redidual = nonlinearstiffness_exp(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, U)

% define function within scope

redidual = -F_forresidual + Ksub_linear * U + Ksub_nonlinear * MaxForce * (exp(U)-1);