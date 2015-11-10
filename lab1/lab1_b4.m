function file = lab1_b4
endfunction

function y = TOL
  y = 10^(-4)
endfunction

function y = N_STEPS
  y = 30
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

function y = f_cyclical_X (x, w, c, a)
  g = c + a * sin(pi*w*x)
  y = x^(g)
endfunction

function yp = fp_cyclical_X (x, w, c, a)
  g = c + a * sin(pi*w*x)
  g_p = a * w * cos(pi*w*x)
  y = x^(g)
  yp = y * (g_p*log(x) + g / x)
endfunction

function y = f_g (X)
  y = 2 * sqrt(X)
endfunction

function y = fp_g (X)
  y = 1 / sqrt(X)
endfunction

function solveLagrangeSystem (_f_X, _fp_X, _title)
  T = 1
  dt = T / N_STEPS

  X = ones(N_STEPS+1, 1)
  l = zeros(N_STEPS, 1)
  alfa_values = zeros(N_STEPS, 1)

  vector_t = linspace(0, T, N_STEPS)

  hold off
  newplot()
  hold on

  i = 0
  while (i++ <= 50)
    n = N_STEPS

    prev_x = X

    l(n) = fp_g(X(n))

    do
      l(n-1) = l(n) * (1 + dt * _fp_X(X(n)) )
    until (--n == 1)

    alfa = @(_n) (l(_n))^(-3/5)

    alfa_values(n) = alfa(n)
    do
      a_n = alfa_values(n+1) = alfa(n+1)
      X(n+1) = X(n) + dt * (_f_X(X(n)) - a_n)
    until (++n == N_STEPS)

    if norm(X - prev_x) < TOL
      break
    end

    plot(vector_t, l, ":m")
    plot(vector_t, alfa_values, ":k")
  end

  plot(vector_t, alfa_values, "-k")
  plot(vector_t, l, "-r")
  title(_title)
endfunction

function doB4()
  solveLagrangeSystem(@f_Xb2, @fp_Xb2, "f(X) from B2")
  print -depsc lab1_b4_1.eps

  solveLagrangeSystem(@f_X2, @fp_X2, "f(X) = x + x^2/10")
  print -depsc lab1_b4_2.eps

  solveLagrangeSystem(@f_X, @fp_X, "f(X) = x")
  print -depsc lab1_b4_3.eps

  w = 10; c=1; a=0.3
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) fp_cyclical_X(x, w, c, a),
		      "f(X): cyclical, w=10, c=1, a=0.3")
  print -depsc lab1_b4_cyclical1.eps

  w = 10; c = 0.1; a = 0.5
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) fp_cyclical_X(x, w, c, a),
		      "f(X): cyclical2, w=10, c=0.1, a=0.5")
  print -depsc lab1_b4_cyclical2.eps

  w = 15; c = 0.5; a = 0.2
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) fp_cyclical_X(x, w, c, a),
		      "f(X): cyclical3, w=15, c=0.5, a=0.2")
  print -depsc lab1_b4_cyclical3.eps

  w = 20
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) (fp_cyclical_X(x, w, c, a)),
		      "f(X): cyclical4, w=20, c=0.5, a=0.2")
  print -depsc lab1_b4_cyclical4.eps

  w = 5; c = 0.1; a = 0.5
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) fp_cyclical_X(x, w, c, a),
		      "f(X): cyclical5, w=5, c=0.1, a=0.5")
  print -depsc lab1_b4_cyclical5.eps

  w = 10; c = -0.3; a = 0.3
  solveLagrangeSystem(@(x) f_cyclical_X(x, w, c, a),
		      @(x) fp_cyclical_X(x, w, c, a),
		      "f(X): cyclical6, w=10, c=-0.3, a=0.3")
  print -depsc lab1_b4_cyclical6.eps
endfunction

doB4()
