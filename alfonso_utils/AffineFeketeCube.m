function affine_fekete_cube = AffineFeketeCube(variables, max_degree, box_scale)
    affine_fekete_cube = FeketeCube(size(variables, 1), max_degree/2, variables);
    affine_fekete_cube = scale_fekete_cube(affine_fekete_cube, box_scale);
end