function cube = sub_vars_fekete_cube(fekete_cube, variables)
    assert(all(size(fekete_cube.mon_basis.variables) == size(variables)));
    cube = fekete_cube;
    cube.polynomials = subs(cube.polynomials, cube.mon_basis.variables, variables);
    cube.mon_basis = sub_vars_mon_basis(cube.mon_basis, variables);

    [P0_to_mon, P_to_mon, mon_to_P0, mon_to_P] = monomial_to_interpolant(cube.P0_large, cube.P_large, cube.polynomials, cube.mon_basis);
    
    cube.P0_to_mon = P0_to_mon;
    cube.P_to_mon = P_to_mon;
    cube.mon_to_P0 = mon_to_P0;
    cube.mon_to_P = mon_to_P;
    
end