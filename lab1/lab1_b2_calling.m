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
    u_v = i * stepsize;   x(i) = u_v
    ls_y(i) = ls_coeffs(1) * u_v + ls_coeffs(2) * u_v^2
    p_y(i) = polyval(polynomial_coeffs, u_v)
  end
  plot(x, ls_y);   hold on;   plot(x, p_y, "--r") #...
endfunction
