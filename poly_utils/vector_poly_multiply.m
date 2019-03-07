function mult_matrix = vector_poly_multiply(poly, in_basis, out_basis)
    %Takes:
    %   poly: msspoly
    %       Written in the variables in `in_basis` and `out_basis`
    %   in_basis: monomial_basis 
    %       as returned by monomial_basis()
    %       This is the basis in which the polynomial on which
    %       `mult_matrix` acts is expressed
    %   out_basis: monomial_basis
    %       as returned by monomial_basis()
    %       This is the basis into which mult_matrix acts
    %       Must be in the same variables as in_basis and
    %       out_basis.degree must be at least in_basis.d + deg(poly)
    %Returns:
    %   mult_matrix: out_basis.num_monomials * in_basis.num_monomials real matrix
    %       mult_matrix * decision_poly_vector_in_in_basis = decision_poly * fixed_poly
    A = sparse(size(out_basis.monomials, 1), size(in_basis.monomials, 1));
    [vars_poly, pows_poly, coeffs_poly] = decomp(poly);
    if size(pows_poly, 2) == 0
        pows_poly = zeros(1, size(out_basis.variables, 1));
        if numel(coeffs_poly) == 0
            coeffs_poly = zeros(1, 1);
        end
    end


    for col_ind=1:size(A, 2)
        [vars_mon, pows_mon, coeffs_mon] = decomp(in_basis.monomials(col_ind));
        for term_ind=1:size(pows_poly, 1)
            if size(pows_mon, 2) == 0
                pows_mon = zeros(1, size(out_basis.variables, 1));
                if numel(coeffs_mon) == 0
                    coeffs_mon = zeros(1, 1);
                end
            end
            ind_poly = monomial_to_index(pows_poly(term_ind,:), vars_poly, out_basis);
            ind_mon = monomial_to_index(pows_mon(1,:), vars_mon, out_basis);
            row_poly = out_basis.power_matrix(ind_poly, :);
            row_mon = out_basis.power_matrix(ind_mon, :);

            [~, index] = ismember(row_poly + row_mon, out_basis.power_matrix, 'rows');
            if index == 0
                error('Product matrix cannot be represented using `out_basis`')
            end
            A(index, col_ind) = coeffs_poly(1, term_ind);
        end
    end
    mult_matrix = A;
end