function doA2a()
  f_slope = @(x) -fProduction(x)   # fProduction(x) - 2*fProduction(x)

  [ t X ] = computeEulerFixpointApproximation(f_slope, 0, 1, 1)
  plot (t, X)
  print -deps lab1_a2_a.eps
endfunction
