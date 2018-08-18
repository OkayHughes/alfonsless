classdef AffineFeketeBasis < FeketeBasis

methods
    function aff_cube = AffineFeketeBasis(variables, max_degree, box_scale)

        name = 'AffineFeketeCube';
            
        aff_cube = aff_cube@FeketeBasis(size(variables, 1), max_degree);
        aff_cube = aff_cube.scale_fekete_cube(box_scale);
        aff_cube = aff_cube.sub_vars_fekete_cube(variables);
        aff_cube.name = 'AffineFeketeCube';

    end

    function cube = scale_fekete_cube(cube, bounds)
        cube.bounds = bounds;

        lb = bounds(:, 1);
        ub = bounds(:, 2);
        scale   = (ub-lb)/2;
        shift   = (lb+ub)/2;
        
        cube.pts = bsxfun(@plus,bsxfun(@times,cube.pts,scale'),shift');
        %P' should equal P_test

        change_of_variables = (cube.variables - shift) .* scale.^(-1);
        cube.polynomials = subs(cube.polynomials, ...
                                       cube.variables, ...
                                       change_of_variables);

        err = dmsubs(cube.polynomials, cube.variables, cube.pts')' - cube.P0_full;
        if abs(mean(mean(err))) > 1E-8
            warning('something may be wrong with polynomial composition\n average subtitution error = %d', mean(mean(err)))
        end
        %scaled_cube.P_full = P_large;
        %scaled_cube.P0_full = P_stupid;

        %prod(scale) = det(Df), where f is the affine change of variables
        % taking cube.pts to scaled_cube.pts
        cube.w = cube.w * prod(scale); 
    end

    function cube = sub_vars_fekete_cube(cube, variables)
        assert(all(size(cube.variables) == size(variables)));
        cube.polynomials = subs(cube.polynomials, cube.variables, variables);
        cube.variables = variables;

    end

    function evaluations = eval_weight_polys(cube, pts)
        if size(pts, 2) ~= size(cube.vars)
            error('cube is %d dimensional, but pts is %d dimensional', size(pts, 2), size(cube.vars));
        end
        lb = cube.bounds(:, 1);
        ub = cube.bounds(:, 2);
        evaluations = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);

    end

end

end