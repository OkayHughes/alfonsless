function A = A_dec_var(pow_mat, coeff_mat, full_vars,  in_basis, out_basis)
    coeff_mat = coeff_mat';
    coeff_vars = in_basis.coeff_vars;

    ind_full_vars = @(variable) var_to_index(variable, full_vars);

    coeff_var_inds = arrayfun(ind_full_vars, coeff_vars);
    coeff_ind_mask = find(coeff_var_inds);
    coeff_var_inds = coeff_var_inds(coeff_ind_mask);
    coeff_stupid = 1:size(coeff_vars, 1);
    coeff_stupid = coeff_stupid(coeff_ind_mask);

    out_var_inds = arrayfun(ind_full_vars, out_basis.variables);
    out_ind_mask = find(out_var_inds);
    out_var_inds = out_var_inds(out_ind_mask);
    out_stupid = 1:size(out_basis.variables, 1);
    out_stupid = out_stupid(out_ind_mask);
    
    pow_mat_coeff_vars = zeros(size(pow_mat, 1), size(coeff_vars, 1));

    pow_mat_out_vars = zeros(size(pow_mat, 1), size(out_basis.variables, 1));

    pow_mat_coeff_vars(:, coeff_stupid) = pow_mat(:, coeff_var_inds);
    pow_mat_out_vars(:, out_stupid) = pow_mat(:, out_var_inds);

    assert(all(sum(pow_mat_coeff_vars, 2) <=1));

    row_inds = sum(pow_mat_coeff_vars, 2) == 1;

    pow_mat_coeff_vars = pow_mat_coeff_vars(row_inds, :);
    pow_mat_out_vars = pow_mat_out_vars(row_inds, :);
    coeff_mat_filt = coeff_mat(row_inds, :);

    A = zeros(out_basis.num_monomials, in_basis.num_monomials);

    pow_mat_out_inds = zeros(size(pow_mat_out_vars, 1), 1);

    for row_ind=1:size(pow_mat_out_vars, 1)
        pow_mat_out_inds(row_ind, 1) = monomial_to_index(pow_mat_out_vars(row_ind, :), ...
                                                         out_basis.variables, out_basis);
    end

    if ~all(pow_mat_out_inds ~= 0)
        error('the basis passed cannot represent the linear operator given by pow_mat');
    end

    for mon_ind=1:in_basis.num_monomials
        mon_present_inds = find(pow_mat_coeff_vars(:, mon_ind));
        out_mon_inds = pow_mat_out_inds(mon_present_inds);

        A(out_mon_inds, mon_ind) = coeff_mat_filt(mon_present_inds, 1);
    end


end
