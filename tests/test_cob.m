bivariate_interp = FeketeCube(2, 6);
bivariate_mon = monomial_basis(2, 6, 'j');

for i=1:100
    coeff_vec = -100 * (rand(1, size(bivariate_mon.monomials, 1)) - 0.5);
    poly = coeff_vec * bivariate_mon.monomials;
    in_interp = dmsubs(poly, bivariate_mon.variables, bivariate_interp.pts')';
    in_interp_from_matrix = bivariate_interp.mon_to_P0 * coeff_vec';
    dist = norm(in_interp_from_matrix - in_interp);                   h
    if dist > 1E-6
        warning('mon_to_P0 is badly behaved')
    end
    %in_interp = in_interp(1:28, :);
    coeff_vec_primed = bivariate_interp.P0_to_mon * in_interp;
    dist = norm(coeff_vec_primed - coeff_vec');
    if dist > 1E-6
       warning('P0_to_mon is badly behaved')
    end
end