function ind = var_to_index(var, vars)

    [~, var_id] = isfree(var);
    [~, vars_ids] = isfree(vars); 
    ind = find(vars_ids == var_id);
    if numel(ind) == 0
        ind = 0;
    end
end