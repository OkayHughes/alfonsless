function affine_fekete_cube = AffineFeketeCube(variables, max_degree, box_scale)
    name = "FeketeCube";
    if check_for_basis(name, size(variables, 1), max_degree/2)
        fekete_cube = read_basis(name, size(variables, 1), max_degree/2);
    else
        sprintf('Creating new %s with n=%d, d=%d', name, size(variables, 1), max_degree)
        fekete_cube = FeketeCube(size(variables, 1), max_degree/2);
        write_basis(fekete_cube);
    end

    affine_fekete_cube = scale_fekete_cube(fekete_cube, box_scale);
    affine_fekete_cube = sub_vars_fekete_cube(affine_fekete_cube, variables);
end