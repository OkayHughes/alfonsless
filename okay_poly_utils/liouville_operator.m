function L_f = liouville_operator(f, in_basis, out_basis)
    v_to_vf = monomial_to_monomial(in_basis, out_basis);
    time_der = vector_derivative(1, in_basis);
    L_f = v_to_vf * time_der;

    for var_ind=1:size(in_basis.variables, 1)-1
        par_der = vector_derivative(1+var_ind, in_basis);
        mult_mat = vector_poly_multiply(f(var_ind), in_basis, out_basis);
        L_f = L_f + mult_mat * par_der;
    end

end