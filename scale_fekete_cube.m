function scaled_cube = scale_fekete_cube(cube, bounds)
    scaled_cube = cube;
    scaled_cube.bounds = bounds;

    lb = bounds(:, 1);
    ub = bounds(:, 2);
    scale   = (ub-lb)/2;
    shift   = (lb+ub)/2;
    scaled_cube.pts = bsxfun(@plus,bsxfun(@times,scaled_cube.pts,scale'),shift');
    %P' should equal P_test
    P_stupid = dmsubs(scaled_cube.polynomials, ...
                      scaled_cube.mon_basis.variables, ...
                      scaled_cube.pts');
    
    [P_large, ~] = qr(P_stupid);
    
    scaled_cube.P_full = P_large;
    scaled_cube.P0_full = P_stupid;
    
    [P0_to_mon, P_to_mon, mon_to_P0, mon_to_P] = monomial_to_interpolant(scaled_cube.P0_full, scaled_cube.P_full, scaled_cube.polynomials, scaled_cube.mon_basis);
    
    intParams.P0_to_mon = P0_to_mon;
    intParams.P_to_mon = P_to_mon;
    intParams.mon_to_P0 = mon_to_P0;
    intParams.mon_to_P = mon_to_P;
end