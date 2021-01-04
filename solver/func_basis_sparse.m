function [result] = func_basis_sparse(x,indices,dim,AA,AT,mode)

if mode == 1
    %reverse sparse basis x->y
    x = reshape(x,dim);
    array = AT(x);
    result = array(indices);
else
    %sparse basis (inverse = conjugate transpose) y_s->x
    expanded_vector = zeros(dim);
    expanded_vector(indices) = x;
    result = AA( expanded_vector );
    result = reshape(result,prod(dim),1);
end

end