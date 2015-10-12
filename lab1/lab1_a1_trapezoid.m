function trap_sum = computeTrapezoidApproximation (f, x, x_end)
  n = 0;   stepsize = (x_end - x) / N_STEPS

  trap_sum = f(x) / 2
  
  while (n++ < N_STEPS)
    x += stepsize
    trap_sum += f(x) * stepsize
  endwhile
  
  trap_sum += f(x_end) / 2
endfunction
