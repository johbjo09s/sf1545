function testFixpointApproximation()
  f_sqrt_2 = @(x) (x + 2 / x) / 2
  computeFixpointApproximation(f_sqrt_2, 3/2)
endfunction

function value = computeFixpointApproximation (f, value)
  iter = 0 ; delta_f = 1

  while (delta_f > TOL) && (iter++ < 20)
    x = value
    value = f(x)
    delta_f = abs(value - x)
  endwhile
endfunction
