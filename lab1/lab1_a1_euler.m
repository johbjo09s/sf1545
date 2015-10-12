function [x_out y_out] = computeEulerApproximation (f, x, x_end, y)
  n = 1;   x_out(1) = x;   y_out(1) = y

  stepsize = (x_end - x) / N_STEPS
  
  while ((x += stepsize) <= x_end)
    y = y + stepsize * f(x)

    n++;
    x_out(n) = x
    y_out(n) = y
  endwhile
endfunction
