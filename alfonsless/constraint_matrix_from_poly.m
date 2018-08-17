function [As, c] = constraint_matrix_from_poly(pol, dec_var_bases, out_basis)
    num_decs = size(dec_var_bases, 2);
    [vars, pows, coeffs] = decomp(pol);
    if nargin < 3
        %Find appropriate out_basis
        dec_var_coeffs = dec_var_bases(1).coeff_vars;
        for dec_ind=2:num_decs
            dec_var_coeffs = [dec_var_coeffs; dec_var_bases(dec_ind).coeff_vars];
        end
        to_degreeify = msubs(pol, dec_var_coeffs, ones(size(dec_var_coeffs)));
        d = msspoly_degree(to_degreeify);

        space_vars = var_set_minus(vars, dec_var_coeffs);

        n = size(space_vars, 1);

        out_basis = monomial_basis(n, d, space_vars);
    end

    As = cell(1, num_decs);
    for dec_ind=1:num_decs
        As{dec_ind} = A_dec_var(pows, coeffs, vars, ...
                      dec_var_bases(1, dec_ind), out_basis);
    end
    const_poly = msubs(pol, dec_var_coeffs, zeros(size(dec_var_coeffs)));
    c = msspoly_to_vector(const_poly, out_basis);

end