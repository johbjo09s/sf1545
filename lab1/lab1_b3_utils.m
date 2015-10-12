function f_matrix = computeF (F, x_values, c_args)
  f_matrix = zeros(1, rows(x_values))
  for i = 1:rows(x_values)
    args = num2cell([x_values(i,:), c_args(1,:)])
    f_matrix(i) = F(args{:})
  end
  f_matrix = transpose(f_matrix)
endfunction

function jacobian = computeJacobian (partial_derivatives, x_values, c_args)
  jacobian = zeros(rows(x_values), columns(partial_derivatives))
  for j = 1:columns(partial_derivatives)
    f_partial = (partial_derivatives{j})
    for i = 1:rows(x_values)
      args = num2cell([x_values(i,:), c_args(1,:)])
      jacobian(i, j) = f_partial(args{:})
    end
  end
endfunction
