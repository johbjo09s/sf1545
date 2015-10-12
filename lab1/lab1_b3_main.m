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
