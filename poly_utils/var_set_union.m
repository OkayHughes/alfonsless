function uni = var_set_union(larg, rarg)
    isnt_in = @(variable) var_to_index(variable, larg) == 0;
    rinds = arrayfun(isnt_in, rarg);
    uni = [larg; rarg(rinds)];
end