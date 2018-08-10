in_basis = monomial_basis(2, 3, 'a');
out_basis = monomial_basis(2, 6, 'a');
vars = in_basis.variables;
const_polys = [vars(1); msspoly(0); msspoly(1); ...
               vars(1) * vars(2); vars(1)^2; ...
               vars(1)*vars(2) +  ...
               2*vars(1)^3 + 40 * vars(2)^3 + 2.2];

tests = [vars(1); msspoly(1); ...
         vars(1) * vars(2); ...
         vars(1)^2 + 2 * vars(2)^2; ...
         vars(1)*vars(2) * 4 + ...
         vars(1) * vars(2)^2 + ...
         1 + vars(1)^3 + vars(2)];

err = 0;
for i=1:size(const_polys, 1)
    mult_matrix = vector_poly_multiply(const_polys(i), in_basis, out_basis);
    for j=1:size(tests, 1)
        vec = mult_matrix * msspoly_to_vector(tests(j), in_basis);
        matrix_result = vec' * out_basis.monomials;
        direct_result = tests(j) * const_polys(i);

        normm = norm(vec - msspoly_to_vector(direct_result, out_basis))
        err = err + normm;
    end
end

if err > 1
    error("Multiplication matrix is totally fucked")
elseif err > 1E-8
    error("Multiplication matrix is not precise")
end