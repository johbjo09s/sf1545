function [x_values y_values] = computeEulerApproximation (f, x, x_end, value)
  iter = 1
  stepsize = (x_end - x) / N_EULER_STEPS

  x_values(iter) = x
  y_values(iter) = value

  while ((x += stepsize) <= x_end)
    iter++;

    value = value + stepsize * f(x)

    x_values(iter) = x
    y_values(iter) = value
  endwhile
endfunction
