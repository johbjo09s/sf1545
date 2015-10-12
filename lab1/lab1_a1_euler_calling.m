function doA1()
  f_slope = @(x) fProduction(x) / 2
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)

  f_slope = @(x) fProduction(x) - 1/2
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)

  f_slope = @(x) fProduction(x) - 0
  [ t X ] = computeEulerApproximation(f_slope, 0, 1, 1)
  plot (t, X)
endfunction
