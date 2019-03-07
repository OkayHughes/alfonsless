function new_basis = sub_vars_mon_basis(mon_basis, new_vars)
    assert(all(size(mon_basis.variables) == size(new_vars)));

    new_basis = mon_basis;

    new_basis.variables = new_vars;
    new_basis.monomials = subs(mon_basis.monomials, mon_basis.variables, new_vars);
    
end