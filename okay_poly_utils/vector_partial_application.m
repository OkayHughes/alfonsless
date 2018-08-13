function application_matrix = vector_partial_application(dimension, value, basis)
    %Takes:
    %   dimension: integer
    %       Dimension along which the returned operator should apply a
    %       polynomial, e.g. apply f(x1, x2)  along dimension x1
    %   value: double
    %       Value to set variable `dimension` to (see below)
    %Returns:
    %   application_matrix: basis.num_monomials square real matrix
    %           Linear operator that computes f(*,..., value, ..., *), and
    %           takes R[x_1, x_2, ..., x_dimension, ..., x_n] -> 
    %           R[x_1, x_2, ..., x_{dimension-1}, x_{dimension+1}, ..., x_n]
    
    vals = zeros(basis.num_monomials, 1);
    row_inds = zeros(basis.num_monomials, 1);
    col_inds = 1:basis.num_monomials;
    
    for mon_ind=1:basis.num_monomials
        row_rep = basis.power_matrix(mon_ind, :);
        scalar = value ^ row_rep(dimension);
        row_rep(dimension) = 0;
        [~, row_index] = ismember(row_rep, basis.power_matrix, 'rows');
        row_inds(mon_ind) = row_index;
        vals(mon_ind) = scalar;
    end
    
    index = find(row_inds ~= 0);
    sanitized_row_inds = row_inds(index);
    sanitized_col_inds = col_inds(index);
    sanitized_vals = vals(index);
    
    application_matrix = sparse(sanitized_row_inds, sanitized_col_inds, sanitized_vals, basis.num_monomials, basis.num_monomials);
    
end