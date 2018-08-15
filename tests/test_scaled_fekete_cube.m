cube = FeketeCube(2, 4);
err = 0;
sc = scale_fekete_cube(cube, [-2, 3; -2, 2]);
for i=1:size(sc.P0_full, 2)
    f = (sc.P0_to_mon * sc.P0_full(:, i))' * sc.mon_basis.monomials;
    tmp_vec = dmsubs(f, sc.mon_basis.variables, sc.pts');
    err = err + norm(tmp_vec' - sc.P0_full(:, i));
end

if err > 1
    error('scaled COB matrix is hot trash')
elseif err > 1E-6
    warning('Scaled COB matrix may be inaccurate')
end