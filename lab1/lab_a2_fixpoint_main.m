function [x_values y_values] = computeEulerFixpointApproximation (f_slope, x, x_end, y)
  n = 1;   x_values(1) = x;   y_values(1) = y

  stepsize = (x_end - x) / N_STEPS
  
  while ((x += stepsize) <= x_end)
    n++;

    f_inner = @(y_next) y + stepsize * f_slope(y_next)
    y = computeFixpointApproximation(f_inner, y)

    x_values(n) = x
    y_values(n) = y
  endwhile
endfunction
