function B3_1 ()
  a_values = U_DATA() (:,1)
  u_values = U_DATA() (:,2)

  U    = @(alfa, a) (8*alfa) / (alfa + 8*a)
  U_a  = @(alfa, a) (8*alfa) * (-8/(alfa + 8*a)^2)

  U_partials = {U_a}

  start_c = [56]

  _compute_U = @(ck) computeF(U, a_values, ck) - u_values
  _compute_J = @(ck) computeJacobian(U_partials, a_values, ck)

  solution = computeGaussNewtonFit(_compute_U, _compute_J, start_c)

  U_fun = @(alfa) U(alfa, solution(1))

  for  i = 1:30
    u_v = i*200;   x(i) = u_v;    y(i) = U_fun(u_v)
  end

  plot(x, y);
endfunction
