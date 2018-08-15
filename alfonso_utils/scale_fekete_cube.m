function scaled_cube = scale_fekete_cube(cube, bounds)
    scaled_cube = cube;
    scaled_cube.bounds = bounds;

    lb = bounds(:, 1);
    ub = bounds(:, 2);
    scale   = (ub-lb)/2;
    shift   = (lb+ub)/2;
    
    scaled_cube.pts = bsxfun(@plus,bsxfun(@times,scaled_cube.pts,scale'),shift');
    %P' should equal P_test

    change_of_variables = (scaled_cube.mon_basis.variables - shift) .* scale.^(-1);
    scaled_cube.polynomials = subs(scaled_cube.polynomials, ...
                                   scaled_cube.mon_basis.variables, ...
                                   change_of_variables);

    err = dmsubs(scaled_cube.polynomials, scaled_cube.mon_basis.variables, scaled_cube.pts')' - scaled_cube.P0_full;
    if norm(err) > 1E-5
        warning('something may be wrong with polynomial composition\n average subtitution error = %d', mean(mean(err)))
    end
    %scaled_cube.P_full = P_large;
    %scaled_cube.P0_full = P_stupid;

    %prod(scale) = det(Df), where f is the affine change of variables
    % taking cube.pts to scaled_cube.pts
    scaled_cube.w = scaled_cube.w * prod(scale); 
end