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
