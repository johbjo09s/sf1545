function doA2b()
#  f_slope = @(x) fProduction(x) - 2*fProduction(x)

  f_slope = @(x) -fProduction(x)
  f_slope_prime = @(x) -1 -(2*x)/5

  [ t X ] = computeEulerNewtonApprox(f_slope, f_slope_prime, 0, 1, 1)
  plot (t, X)
  print -deps lab1_a2_b.eps
endfunction
