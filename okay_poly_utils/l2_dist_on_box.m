function dst = l2_dist_on_box(p1, p2, lb, ub, vars)
    integrand = (p1-p2)^2;
    [vars, ~, ~] = decomp(integrand);
    for i=1:size(vars, 1)
        integrand = subs(integral(integrand,vars(i)),vars(i),ub(i)) ...
            - subs(integral(integrand,vars(i)),vars(i),lb(i));
    end
    dst = dmsubs(integrand, vars, zeros(size(vars)));
end