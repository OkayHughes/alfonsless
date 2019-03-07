function poly_vector = msspoly_to_vector(poly, basis)
    %converts a mss_poly expressed in terms of the variables in
    %basis.variables into a vector corresponding to basis.monomials,
    %i.e, poly_vector' * basis.monomials = poly.
    poly_vector = zeros(basis.num_monomials, 1);
    [vars, pows, coeffs] = decomp(poly);
    
    %handle edge case for constant polynomials (ugh)
    if size(pows, 2) == 0
        idx = monomial_to_index(msspoly(1), vars, basis);
        if numel(coeffs) == 0
            poly_vector(idx) = 0;
        else
            poly_vector(idx) = coeffs(1, 1);
        end
        return
    end
    
    for term_ind=1:size(coeffs, 2)
       term_degree = full(sum(pows(term_ind)));
       if term_degree > basis.d
           error('Error:\nPolynomial has term of degree %d, but basis has max degree %d', term_degree, basis.d)
       end
       idx = monomial_to_index(pows(term_ind, :), vars, basis);
       if idx == 0
            error('Error:\nPolynomial has term containing a monomial not in basis.monomials');
       end
       poly_vector(idx) = coeffs(1, term_ind);
    end
end