function res = var_set_eq(larg, rarg)
    res = var_set_subset(larg, rarg) & var_set_subset(larg, rarg);
end