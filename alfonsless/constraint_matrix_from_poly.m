%NOTE THIS FUNCTION CONTAINS A HACK
function [As, c, out_basis] = constraint_matrix_from_poly(pol, dec_var_bases, variables, degree)
    num_decs = size(dec_var_bases, 1);
    [vars, pows, coeffs] = decomp(pol);

    %Find appropriate out_basis
    dec_var_coeffs = dec_var_bases(1).coeff_vars;
    for dec_ind=2:num_decs
        dec_var_coeffs = [dec_var_coeffs; dec_var_bases(dec_ind).coeff_vars];
    end

    space_vars = variables;

    %handling the zero-d case. THIS IS A HACK
    if size(space_vars, 1) == 0
        space_vars = dec_var_bases(1).mon_basis.variables(1);
    end

    out_basis = MonomialBasis(space_vars, degree);

    As = cell(1, num_decs);
    for dec_ind=1:num_decs
        As{dec_ind} = A_dec_var(pows, coeffs, vars, ...
                      dec_var_bases(dec_ind), out_basis);
    end
    const_poly = subs(pol, dec_var_coeffs, zeros(size(dec_var_coeffs)));
    c = msspoly_to_vector(const_poly, out_basis);

end