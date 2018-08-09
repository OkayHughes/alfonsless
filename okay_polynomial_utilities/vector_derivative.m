function derivative_operator = vector_derivative(partial_i, basis)
    %Takes:
    %   partial_i: integer
    %       which variable number (indexed by the variables in `basis`)
    %       to differentiate with respect to.
    %   basis: struct
    %       struct returned from `monomial_basis`
    %Returns:
    %   derivative_operator: nchoosek(basis.n + basis.d, basis.d)^2 sparse_matrix 
    %       If a polynomial p(x_1, ..., x_n) is represented by a vector q in
    %       the basis defined by `basis`, then \partial_i(p) =
    %       derivative_operator * q
    
    row_indices = zeros(basis.num_monomials, 1);
    col_indices = 1:basis.num_monomials;
    vals = zeros(basis.num_monomials, 1);
    for term_ind=1:basis.num_monomials
        monomial_coefficients = basis.power_matrix(term_ind, :);
        relevant_power = monomial_coefficients(partial_i);
        if relevant_power-1 < 0
            continue;
        end
        monomial_coefficients(partial_i) = relevant_power-1;
        %note: this has a bunch of dumb linear searches, might be able to sort
        %matrix in a good way
        
        %TODO: some kind of stupid ass indexing bug here
        [~, result_polynomial_ind] =  ismember(monomial_coefficients, basis.power_matrix, 'rows');
        %row_indices
        row_indices(term_ind) = result_polynomial_ind;
        vals(term_ind) = relevant_power;
    end
    
    index = find(row_indices ~= 0);
    sanitized_row_indices = row_indices(index);
    sanitized_col_indices = col_indices(index);
    sanitized_vals = vals(index);
    derivative_operator = sparse(sanitized_row_indices, sanitized_col_indices, sanitized_vals, basis.num_monomials, basis.num_monomials);
end