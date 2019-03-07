function integ = def_int_on_box(p, lb, ub, vars)

    integrand = p;
    [vars, ~, ~] = decomp(integrand);
    for i=1:size(vars, 1)
        integrand = subs(integral(integrand,vars(i)),vars(i),ub(i)) ...
            - subs(integral(integrand,vars(i)),vars(i),lb(i));
    end
    integ = dmsubs(integrand, vars, zeros(size(vars)));
end