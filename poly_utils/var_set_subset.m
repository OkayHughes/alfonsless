function res = var_set_subset(larg, rarg)
    is_in = @(variable) var_to_index(variable, rarg) ~= 0;
    res = all(arrayfun(is_in, larg));

end