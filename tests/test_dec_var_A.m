
v = msspoly('x', 3);
in_basis = monomial_basis(2, 4, v(1:2, 1), msspoly('@a', nchoosek(2 + 4, 2)));
out_basis = monomial_basis(3, 6, v, msspoly('@b', nchoosek(3 + 4, 3)));
f = v(1)^2 + v(3)^2;
A = vector_poly_multiply(f, in_basis, out_basis);
pol = transpose(in_basis.coeff_vars) * in_basis.monomials; 
pr = f * pol;

[vars, pows, coeffs] = decomp(pr);
A_new = A_dec_var(pows, coeffs, vars, in_basis, out_basis);

normm = norm(A - A_new)

v = msspoly('x', 2);
q = msspoly('y', 2);
in_basis = monomial_basis(2, 4, v, msspoly('@a', nchoosek(2 + 4, 2)));
out_basis = monomial_basis(4, 6, [v;q], msspoly('@b', nchoosek(4 + 6, 4)));
f = v(1)^2 + q(2)^2;
A = vector_poly_multiply(f, in_basis, out_basis);
pol = transpose(in_basis.coeff_vars) * in_basis.monomials; 
pr = f * pol;

[vars, pows, coeffs] = decomp(pr);
A_new = A_dec_var(pows, coeffs, vars, in_basis, out_basis);

normm = norm(A - A_new)

v = msspoly('x', 3);
q = msspoly('y', 2);
in_basis = monomial_basis(5, 6, [v;q], msspoly('@a', nchoosek(5 + 6, 5)));
out_basis = monomial_basis(2, 6, [v(1);q(1)], msspoly('@b', nchoosek(2 + 6, 2)));

l2s = monomial_to_monomial(in_basis, out_basis);
partial_v2 = vector_partial_application(2, 0.5, in_basis);
partial_v3 = vector_partial_application(3, 0, in_basis);
partial_q2 = vector_partial_application(5, -1, in_basis);
A = l2s * partial_v2 * partial_v3 * partial_q2;

pol = transpose(in_basis.coeff_vars) * in_basis.monomials; 
app = msubs(pol, v(2), 0.5);
app = msubs(app, v(3), 0);
app = msubs(app, q(2), -1);

[vars, pows, coeffs] = decomp(app);

A_new = A_dec_var(pows, coeffs, vars, in_basis, out_basis);

normm = norm(A - A_new)


v = msspoly('x', 3);
in_basis = monomial_basis(3, 6, [v], msspoly('@a', nchoosek(3 + 6, 3)));
arb_basis = monomial_basis(3, 6, [v], msspoly('@c', nchoosek(3 + 6, 3)));
out_basis = monomial_basis(3, 6, [v], msspoly('@b', nchoosek(3 + 6, 3)));
der_1 = vector_derivative(1, in_basis);
der_2 = vector_derivative(2, in_basis);
A =  der_1 * der_2;

pol = transpose(in_basis.coeff_vars) * in_basis.monomials; 
app = diff(pol, v(1));
app = diff(app, v(2));

pol2 = transpose(arb_basis.coeff_vars) * arb_basis.monomials; 

app = app + diff(pol2, v(2));

[vars, pows, coeffs] = decomp(app);

A_new = A_dec_var(pows, coeffs, vars, in_basis, out_basis);

normm = norm(A - A_new)






