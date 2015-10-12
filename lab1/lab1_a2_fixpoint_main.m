function y = computeNewtonApproximation (f, fprime, start)
  f_newton = @(x) x - (f(x) / fprime(x))
  y = computeFixpointApproximation(f_newton, start)
endfunction
