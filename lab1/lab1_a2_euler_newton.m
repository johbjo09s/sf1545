function [x_out y_out] = computeEulerNewtonApprox (f, fprime, x, x_end, y)
  n = 1;   x_out(1) = x;   y_out(1) = y

  stepsize = (x_end - x) / N_STEPS

  while ((x += stepsize) <= x_end)
    f_euler       = @(e_x) y + stepsize * f(e_x)      - e_x 
    f_euler_prime = @(e_x)     stepsize * fprime(e_x) - 1

    y = computeNewtonApproximation(f_euler, f_euler_prime, y)

    n++;
    x_out(n) = x
    y_out(n) = y
  endwhile
endfunction
