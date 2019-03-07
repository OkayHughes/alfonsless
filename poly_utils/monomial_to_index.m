function index = monomial_to_index(power_matrix_row, variables, basis)
    if size(variables, 1) == 0
        row_rep = zeros(basis.n, 1);
        [~, index] = ismember(row_rep', basis.power_matrix, 'rows');
        return;
    end
    [~, ref_var_ids] = isfree(basis.variables);
    [~, inp_var_ids] = isfree(variables);
    row_rep = zeros(basis.n, 1);
    for var_ind=1:size(power_matrix_row, 2)
        row_ind = find(ref_var_ids == inp_var_ids(var_ind));
        if numel(row_ind) == 0
            index = 0;
            return
        end
        row_rep(find(ref_var_ids == inp_var_ids(var_ind))) = power_matrix_row(var_ind);
    end
    [~, index] = ismember(row_rep', basis.power_matrix, 'rows');
end