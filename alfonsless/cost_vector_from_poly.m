function bs = cost_vector_from_poly(pol, dec_var_polys)
    num_decs = size(dec_var_bases, 2);
    dec_var_coeffs = dec_var_bases(1).coeff_vars;
    for i=2:num_decs
        dec_var_coeffs = [dec_var_coeffs; dec_var_bases(i).coeff_vars];
    end

    [vars, pows, coeffs] = decomp(pol);
    b = zeros(size(dec_var_coeffs));
    for i=1:size(pows, 1)
        assert(sum(pows(i, :), 2) == 1)
        ind = find(pows(i, :));
        v = vars(ind);
        ind = var_to_index(v, dec_var_coeffs);
        b(ind) = coeffs(1, i);
    end

    bs = {};
    tally = 1;
    for i=1:num_decs
        b_i = b(tally:tally - 1 + dec_var_bases(i).mon_basis.num_monomials);
        bs = {bs; b_i};
        tally = tally + dec_var_bases(i).mon_basis.num_monomials;
    end
end