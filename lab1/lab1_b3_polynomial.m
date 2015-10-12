function p = computePolyfit (f_values)
  x = transpose(f_values(:,1))
  y = transpose(f_values(:,2))
  p = polyfit(x, y, rows(f_values))
endfunction

function B3_polynomial ()
  polynomial_coeffs = computePolyfit(U_DATA)

  for  i = 1:30
    u_v = i*200;   x(i) = u_v;   y(i) = polyval(polynomial, u_v)
  end
				   
  plot(x, y)
  hold on
  plot(U_DATA() (:,1), U_DATA() (:,2), "or", "markersize", 5)
  print -deps lab1_b3_polynomial.eps
  hold off
endfunction
