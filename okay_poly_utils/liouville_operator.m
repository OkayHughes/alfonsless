function L_f = liouville_operator(f, time_index, first_sp_ind, in_basis, out_basis)
    v_to_vf = monomial_to_monomial(in_basis, out_basis);
    time_der = vector_derivative(time_index, in_basis);
    L_f = v_to_vf * time_der;
    ct = 1;
    for var_ind=first_sp_ind:(first_sp_ind + size(f, 1) - 1)
        par_der = vector_derivative(var_ind, in_basis);
        mult_mat = vector_poly_multiply(f(ct), in_basis, out_basis);
        L_f = L_f + mult_mat * par_der;
        ct = ct + 1;
    end

end