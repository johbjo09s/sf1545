function file = lab1_b2
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

  print -depsc lab1_b2.eps
  hold off
endfunction

B2()
