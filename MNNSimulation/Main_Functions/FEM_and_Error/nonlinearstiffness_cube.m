function redidual = nonlinearstiffness_cube(F_forresidual, Ksub_linear, Ksub_nonlinear, MaxForce, MaxU, U)

% define function within scope

redidual = -F_forresidual + Ksub_linear * U + Ksub_nonlinear * MaxForce * ((U./1).^3);