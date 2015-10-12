  % Integrering
  f_alfa = @(x) 1/2
  U      = @(x) - 1 / f_alfa(x)
  int_U  = @(x) f_U(x) / (fProduction(x) - f_alfa(x))
  int_U_approx = computeTrapezoidApproximation(int_U, 0, 1)  
