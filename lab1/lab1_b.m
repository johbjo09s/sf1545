
function b = lab1_b
endfunction

function y = TOL
  y = 10^(-4)
endfunction

function y = N_STEPS
  y = 10
endfunction

function y = fProduction (x)
  y = x + ((x^2) / 5)
endfunction

function [solution, residual] = computeLeastSquaresFit (varargin)
  assert(nargin > 2)

  f_values = (varargin{1})

  num_terms = nargin - 1

  y = zeros(rows(f_values), 1)
  A = zeros(rows(f_values), num_terms)

  for this_row = 1:rows(f_values)
    y(this_row, 1) = f_values(this_row, 2)
    for this_column = 2:(nargin)
      f_term = (varargin{this_column})
      f_arg = f_values(this_row, 1)
      A(this_row, this_column-1) = f_term(f_arg)
    end
  end

  solution = A \ y
  residual = y - A*solution
endfunction

function testLeastSquaresFit ()
  f_values = [40   , 55.3  ;
	      45   , 71.9  ;
	      50   , 92.5  ;
	      55   , 118   ;
	      60   , 149.4 ]

  [ coeffs, residual ] = computeLeastSquaresFit(f_values, @(x) 1, @(x) x)

  stepsize = (max(f_values(:,1) - min(f_values(:,1)))) / 50

  for i = 1:50
    u_v = i * stepsize
    x(i) = u_v
    ls_y(i) = coeffs(1) * 1 + coeffs(2) * u_v
    #p_y(i) = polyval(polynomial_coeffs, u_v)
  end

  plot(x, ls_y)
  hold on
  plot(f_values(:,1), f_values(:,2), "or", "markersize", 5)
  # plot(f_values(:,1), residual, "or", "markersize", 5)
  pause
endfunction

#testLeastSquaresFit()
#exit()

function solution = computeGaussNewtonFit (_compute_F, _compute_J, ck)
  iter = 0
  delta_f = 1

  while (delta_f > TOL) && (iter++ < 30)
    c_prev = ck
    
    f = _compute_F(ck)
    jacobian = _compute_J(ck)
 
    delta_c = jacobian \ f

    ck = ck - transpose(delta_c)

    delta_f = norm(ck - c_prev)
  endwhile

  solution = ck
endfunction

function f_matrix = computeF (F, x_values, c_args)
  f_matrix = zeros(1, rows(x_values))
  for i = 1:rows(x_values)
    args = num2cell([x_values(i,:), c_args(1,:)])
    f_matrix(i) = F(args{:})
  end
  f_matrix = transpose(f_matrix)
endfunction

function jacobian = computeJacobian (partial_derivatives, x_values, c_args)
  jacobian = zeros(rows(x_values), columns(partial_derivatives))
  for j = 1:columns(partial_derivatives)
    f_partial = (partial_derivatives{j})
    for i = 1:rows(x_values)
      args = num2cell([x_values(i,:), c_args(1,:)])
      jacobian(i, j) = f_partial(args{:})
    end
  end
endfunction

function testGaussNewtonFit1 ()
  f_values = [0.0 , 7.4 ;
	      0.1 , 6.6 ;
	      0.2 , 6.0 ;
	      0.3 , 5.6 ;
	      0.4 , 5.2 ;
	      0.5 , 4.9 ;
	      0.6 , 4.7 ;
	      0.7 , 4.4 ;
	      0.8 , 4.2 ;
	      0.9 , 4.1 ;
	      1.0 , 4.0 ]

  x_values = f_values(:,1)
  y_values = f_values(:,2)

  F    = @(x, c1, c2, c3) c1 + c2*e^(-c3*x)
  F_c1 = @(x, c1, c2, c3) 1
  F_c2 = @(x, c1, c2, c3)         e^(-c3*x)
  F_c3 = @(x, c1, c2, c3)  - x*c2*e^(-c3*x)

  F_partials = {F_c1, F_c2, F_c3}

  _compute_F = @(ck) computeF(F, x_values, ck) - y_values
  _compute_J = @(ck) computeJacobian(F_partials, x_values, ck)

  start_c = [3.5, 3.9, 2]

  solution = computeGaussNewtonFit(_compute_F, _compute_J, start_c)
  residual = y_values - F([x_values, solution]{:})
endfunction

#testGaussNewtonFit1()

function testGaussNewtonFit2 ()
  f_values = [0.2 , 3.16 ;
	      0.3 , 2.38 ;
	      0.4 , 1.75 ;
	      0.5 , 1.34 ;
	      0.6 , 1.00 ]

  x_values = f_values(:,1)
  y_values = f_values(:,2)

  F   = @(x, a, b)    a*e^(-b*x)
  F_a = @(x, a, b)      e^(-b*x)
  F_b = @(x, a, b) -x*a*e^(-b*x)

  F_partials = {F_a, F_b}

  _compute_F = @(ck) computeF(F, x_values, ck) - y_values
  _compute_J = @(ck) computeJacobian(F_partials, x_values, ck)

  start_c = [5.618, 2.88]

  solution = computeGaussNewtonFit(_compute_F, _compute_J, start_c)

  # expect = [5.6310 , 2.8872]
endfunction

#testGaussNewtonFit2()

function testGaussNewtonFit3 ()
  f_values = [0.5 , 0.3 ;
	      0.8 , 0.3 ;
	      1.0 , 0.5 ;
	      1.2 , 0.9 ;
	      1.5 , 1.4 ;
	      1.8 , 1.1 ;
	      2.0 , 0.5 ;
	      2.4 , 0.3 ]

  x_values = f_values(:,1)
  y_values = f_values(:,2)

  F    = @(x, a, b, w, x0) a        + b * sin(w*(x - x0))
  F_a  = @(x, a, b, w, x0) 1
  F_b  = @(x, a, b, w, x0)                sin(w*(x - x0))
  F_w  = @(x, a, b, w, x0) (x - x0) * b * cos(w*(x - x0))
  F_x0 = @(x, a, b, w, x0)       -w * b * cos(w*(x - x0))

  F_partials = {F_a, F_b, F_w, F_x0}

  _compute_F = @(ck) computeF(F, x_values, ck) - y_values
  _compute_J = @(ck) computeJacobian(F_partials, x_values, ck)

  start_c = [0.7, 0.7, pi, 1.2]

  solution = computeGaussNewtonFit(_compute_F, _compute_J, start_c)

  # expect = [0.7761, 0.5860, 3.9225, 1.1092]
endfunction

#testGaussNewtonFit3()

function p = computePolyfit (f_values)
  x = transpose(f_values(:,1))
  y = transpose(f_values(:,2))
  p = polyfit(x, y, rows(f_values))
endfunction

function B2 ()
  f_values = [0    , 0    ;
	      0.5  , 0.52 ;
	      1    , 1.09 ;
	      1.5  , 1.75 ;
	      2    , 2.45 ;
	      2.99 , 3.5  ;
	      3    , 4.0]  

  [ls_coeffs, ls_residual] = computeLeastSquaresFit(f_values, @(x) x, @(x) x^2)
  
  polynomial_coeffs = computePolyfit(f_values)

  stepsize = (max(f_values(:,1) - min(f_values(:,1)))) / 50

  for i = 1:50
    u_v = i * stepsize
    x(i) = u_v
    ls_y(i) = ls_coeffs(1) * u_v + ls_coeffs(2) * u_v^2
    p_y(i) = polyval(polynomial_coeffs, u_v)
  end

  plot(x, ls_y)
  hold on
  plot(x, p_y, "--r")
  plot(f_values(:,1), f_values(:,2), "or", "markersize", 5)
  plot(f_values(:,1), ls_residual, "ox", "markersize", 5)
  title(["ax + bx^2, a=" num2str(ls_coeffs(1)) ", b=" num2str(ls_coeffs(2))  ])

  print -deps lab1_b2.eps
  hold off
endfunction

B2()

function v = U_DATA ()
  v = [150  , 2 ;
       200  , 3 ;
       300  , 4 ;
       500  , 5 ;
       1000 , 6 ;
       2000 , 7 ]
#      10*^9 , 8 ]
endfunction

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
    u_v = i*200;   x(i) = u_v;   y(i) = U_fun(u_v)
  end
  plot(x, y);
  title(["8*alfa / (alfa + 8a), a=" num2str(solution(1))])
  hold on
  plot(a_values, u_values, "or", "markersize", 5)
  print -deps lab1_b3_1.eps
  hold off
endfunction

B3_1()

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
  title(["8 - a*(alfa^b), a=" num2str(solution(1)) ", b=" num2str(solution(2))])
  hold on
  plot(a_values, u_values, "or", "markersize", 5)
  print -deps lab1_b3_2.eps
  hold off
endfunction

B3_2()

function B3_polynomial ()
  polynomial = computePolyfit(U_DATA)
  
  for  i = 1:30
    u_v = i*200;   x(i) = u_v;   y(i) = polyval(polynomial, u_v)
  end
  plot(x, y)
  hold on
  plot(U_DATA() (:,1), U_DATA() (:,2), "or", "markersize", 5)
  print -deps lab1_b3_polynomial.eps
  hold off
endfunction

B3_polynomial()


