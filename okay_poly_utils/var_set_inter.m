function inter = var_set_inter(larg, rarg)
    is_in = @(variable) var_to_index(variable, rarg) ~= 0;
    var_inds = arrayfun(is_in, larg);
    inter = larg(var_inds);

end