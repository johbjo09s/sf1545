
function l = lab1_b4
endfunction

function y = TOL
  y = 10^(-4)
endfunction

function y = N_STEPS
  y = 50
endfunction

function y = f_Xb2 (x)
  a = 1.0709
  b = 0.061849
  y = a*x + b*(x^2)
endfunction

function y = fp_Xb2 (x)
  a = 1.0709
  b = 0.061849
  y = a + 2*b*(x)
endfunction

function y = f_X2 (x)
  y = x + ((x^2) / 10)
endfunction

function y = fp_X2 (x)
  y = 1 + (2*x / 10)
endfunction

function y = f_X (x)
  y = x
endfunction

function y = fp_X (x)
  y = 1
endfunction

function y = f_U (alfa)
  y = -(3/2)*alfa^(-2/3)
endfunction

function y = f_g (X)
  y = 2*sqrt(X)
endfunction

function y = fp_g (X)
  y = 1 / sqrt(X)
endfunction

function solveLagrangeSystem (my_f_X, my_fp_X, my_title)
  T = 1

  dt = T / N_STEPS

  X = zeros(N_STEPS, 1)
  l = zeros(N_STEPS, 1)
  alfa_values = zeros(N_STEPS+1, 1)

  n = 0
  while (n++ <= N_STEPS)
    X(n) = 1
  end

  vector_t = linspace(0, T, N_STEPS+1)
  vector_t2 = linspace(0, T, N_STEPS)
  
  i = 0
  while (i++ <= 20)
    n = N_STEPS

    prev_x = X

    l(n+1) = fp_g(X(n))
    l(n) = fp_g(X(n))

    do 
      l_n = l(n)
      l(n-1) = l_n + dt * my_fp_X(X(n)) * l_n
    until (n-- == 2)

    n = 1
    
    alfa_n = @(tmp_n) (l(tmp_n+1))^(-3/5)
#    utility_n = @(tmp_n) - (3/2) * alfa_n(tmp_n)^(2/3)
    
    do
      alfa_values(n+1) = alfa_n(n)
#      utility_values(n+1) = utility_n(n)
      X(n+1) = X(n) + dt * (my_f_X(X(n)) - alfa_n(n))
    until (n++ == N_STEPS)

    dX = X - prev_x

    if norm(dX) < TOL
      i = 100
    end
   # plot(vector_t, l, "--r")
  end

  plot(vector_t, X)
  hold on
  plot(vector_t, alfa_values, "--")
  title(my_title)
endfunction

function doB4()
  solveLagrangeSystem(@f_Xb2, @fp_Xb2, "f(X) from B2")
  print -deps lab1_b4_1.eps
  hold off

  solveLagrangeSystem(@f_X2, @fp_X2, "f(X) = x + x^2/10")
  print -deps lab1_b4_2.eps
  hold off

  solveLagrangeSystem(@f_X, @fp_X, "f(X) = x")
  print -deps lab1_b4_3.eps
  hold off
endfunction

doB4()
