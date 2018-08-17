num_dec_polys = 2;

t = msspoly('t', 1);
x = msspoly('x', 2);

q = [t;x];

in_basis_Z = monomial_basis(3, 6, q, msspoly('@a', nchoosek(3 + 6, 3)));
in_basis_Y = monomial_basis(2, 6, x, msspoly('@b', nchoosek(2 + 6, 2)));

c = 8 * x(1)^2 + 1.5 * x(2) ^ 5 + 9;

part_app_t = vector_partial_application(1, 2, in_basis_Z);
Z_to_Y = monomial_to_monomial(in_basis_Z, in_basis_Y);

A_old = {Z_to_Y * part_app_t, eye(in_basis_Y.num_monomials)};
c_old = msspoly_to_vector(c, in_basis_Y);

v = transpose(in_basis_Z.coeff_vars) * in_basis_Z.monomials; 
vp = msubs(v, t, 2);

w = transpose(in_basis_Y.coeff_vars) * in_basis_Y.monomials; 

pol = vp + w + c;

dec_var_bases = [in_basis_Z, in_basis_Y];

[As, c] = constraint_matrix_from_poly(pol, dec_var_bases);

for i=1:num_dec_polys
    normm_A = norm(As{i} - A_old{i})
end

normm_c = norm(c - c_old)
