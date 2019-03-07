function L_g = liouville_operator_g(g, first_sp_ind, in_basis, out_basis)
    L_g = zeros(out_basis.num_monomials, in_basis.num_monomials);
    ct = 1;
    for var_ind=first_sp_ind:(first_sp_ind + size(g, 1) - 1)
        par_der = vector_derivative(var_ind, in_basis);
        mult_mat = vector_poly_multiply(g(ct), in_basis, out_basis);
        L_g = L_g + mult_mat * par_der;
        ct = ct + 1;
    end

end