function power_mat = monomials_to_matrix(monoms, variables)
    %Takes:
    %   monomials: n x 1 vector of monomials 
    %              like what is returned from monomials() in msspoly
    %   variables: k x 1 vector of free variables
    %              like what is returned from msspoly('x', k)
    %Returns:
    %   power_mat: n x k real matrix A
    %              A_(i,j) is the power of the jth variable in the ith 
    %              monomial in `monomials`
    [~, var_ids] = isfree(variables);
    power_mat = zeros(size(monoms, 1), size(variables, 1));
    for monom_ind = 1:size(monoms, 1)
        [mon_vars, mon_pows, ~] = decomp(monoms(monom_ind));
        [~, mon_var_ids] = isfree(mon_vars);
        for var_ind=1:size(mon_pows, 2)
            %TODO: make sure constant monomial works correctly
            power_mat(monom_ind, find(var_ids == mon_var_ids(var_ind))) = mon_pows(1, var_ind);
        end
    end
end