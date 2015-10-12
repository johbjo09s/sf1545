function B3_2 ()
  a_values = U_DATA() (:,1)
  u_values = U_DATA() (:,2)

  U    = @(alfa, a, b) 8 - a*(alfa^b)
  U_a  = @(alfa, a, b) - alfa^b
  U_b  = @(alfa, a, b) - a * (alfa^b) * log(alfa)

  U_partials = {U_a, U_b}

  _compute_U = @(ck) computeF(U, a_values, ck) - u_values
  _compute_J = @(ck) computeJacobian(U_partials, a_values, ck)

  start_c = [100, -1/2]   % Gissningar, experimentering. Känslig

  solution = computeGaussNewtonFit(_compute_U, _compute_J, start_c)

  U_fun = @(alfa) U(alfa, solution(1), solution(2))

  for  i = 1:30
    u_v = i*200;   x(i) = u_v;   y(i) = U_fun(u_v)
  end
  plot(x, y)
endfunction
