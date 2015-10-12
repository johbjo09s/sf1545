
function a = lab1_a
endfunction

function y = TOL
  y = 10^(-4)
endfunction

function steps = N_STEPS
  steps = 10
endfunction

function y = fProduction (x)
  y = x + ((x^2) / 5)
endfunction

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

function trap_sum = computeTrapezoidApproximation (f, x, x_end)
  n = 0;   stepsize = (x_end - x) / N_STEPS

  trap_sum = f(x) / 2
  
  while (n++ < N_STEPS)
    x += stepsize
    trap_sum += f(x) * stepsize
  endwhile
  
  trap_sum += f(x_end) / 2
endfunction

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

function testNewtonApproximation()
  f = @(x) 2 - x^2
  fprime = @(x) -2*x
  computeNewtonApproximation(f, fprime, 1/6)
endfunction

function y = computeNewtonApproximation (f, fprime, start)
  f_newton = @(x) x - (f(x) / fprime(x))
  y = computeFixpointApproximation(f_newton, start)
endfunction

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

function doA1()
  f_slope = @(x) fProduction(x) / 2
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)
  title("F(X)/2")
  print -deps lab1_a1_i.eps

  f_slope = @(x) fProduction(x) - 1/2
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)
  title("F(X) - 1/2")
  print -deps lab1_a1_ii.eps

  f_slope = @(x) fProduction(x) - 0
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)
  title("F(X) - 0")
  print -deps lab1_a1_iii.eps

  % Integrering
  f_alfa = @(x) 1/2
  f_U    = @(x) - 1 / f_alfa(x)
  int_U  = @(x) f_U(x) / (fProduction(x) - f_alfa(x))
  int_U_approx = computeTrapezoidApproximation(int_U, 0, 1)  
endfunction

function doA2a()
#  f_slope = @(x) fProduction(x) - 2*fProduction(x)
  f_slope = @(x) -fProduction(x)

  [ t X ] = computeEulerFixpointApproximation(f_slope, 0, 1, 1)
  plot (t, X)
  title("-F(X)")
  print -deps lab1_a2_a.eps
endfunction

function doA2b()
#  f_slope = @(x) fProduction(x) - 2*fProduction(x)

  f_slope = @(x) -fProduction(x)
  f_slope_prime = @(x) -1 -(2*x)/5

  [ t X ] = computeEulerNewtonApprox(f_slope, f_slope_prime, 0, 1, 1)
  plot(t, X)
  title("-F(X)")
  print -deps lab1_a2_b.eps
endfunction

function xdot = f_testOde(x, t)
  xdot(1) = -fProduction(t)
endfunction

function TestOde()
  tRange = linspace(0, 1, N_STEPS)

  test_x = lsode(@f_testOde, 1, tRange)

  plot (tRange, test_x)
  print -deps lab1_TestOde.eps
endfunction

function [t_values x_values] = A2_ComputeSolution()
  x_solution = @(a) 5 / (6*exp(a) - 1)
  t = 0
  x_values(1) = 1
  t_values(1) = 0
  iter = 1

  while (t += 0.1) < 1
    x_values(iter) = x_solution(t)
    t_values(iter) = t
    iter++;
  endwhile

  plot (t_values, x_values)
  print -deps lab1_a2_analytic.eps
endfunction

doA1()
doA2a()
doA2b()
