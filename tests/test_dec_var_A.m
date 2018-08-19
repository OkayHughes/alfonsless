
v = msspoly('x', 3);
in_basis = MonomialBasis(v(1:2, 1), 4);

out_basis = MonomialBasis(v, 6);
f = v(1)^2 + v(3)^2;
A = vector_poly_multiply(f, in_basis, out_basis);
pol_obj = Polynomial(in_basis);
pol = pol_obj.to_msspoly();
pr = f * pol;

[vars, pows, coeffs] = decomp(pr);
A_new = A_dec_var(pows, coeffs, vars, pol_obj, out_basis);

normm = norm(A - A_new)

v = msspoly('x', 2);
q = msspoly('y', 2);
in_basis = MonomialBasis(v, 4);
out_basis = MonomialBasis([v;q], 6);
f = v(1)^2 + q(2)^2;
A = vector_poly_multiply(f, in_basis, out_basis);
pol_obj = Polynomial(in_basis);
pol = pol_obj.to_msspoly(); 
pr = f * pol;

[vars, pows, coeffs] = decomp(pr);
A_new = A_dec_var(pows, coeffs, vars, pol_obj, out_basis);

normm = norm(A - A_new)

v = msspoly('x', 3);
q = msspoly('y', 2);
in_basis = MonomialBasis([v;q], 6);
out_basis = MonomialBasis([v(1);q(1)], 6);

l2s = monomial_to_monomial(in_basis, out_basis);
partial_v2 = vector_partial_application(2, 0.5, in_basis);
partial_v3 = vector_partial_application(3, 0, in_basis);
partial_q2 = vector_partial_application(5, -1, in_basis);
A = l2s * partial_v2 * partial_v3 * partial_q2;

pol_obj = Polynomial(in_basis);
pol = pol_obj.to_msspoly();; 
app = msubs(pol, v(2), 0.5);
app = msubs(app, v(3), 0);
app = msubs(app, q(2), -1);

[vars, pows, coeffs] = decomp(app);

A_new = A_dec_var(pows, coeffs, vars, pol_obj, out_basis);

normm = norm(A - A_new)


v = msspoly('x', 3);
in_basis = MonomialBasis([v], 6);
arb_basis = MonomialBasis([v], 6);
out_basis = MonomialBasis([v], 6);
der_1 = vector_derivative(1, in_basis);
der_2 = vector_derivative(2, in_basis);
A =  der_1 * der_2;

pol_obj = Polynomial(in_basis);
pol = pol_obj.to_msspoly();
app = diff(pol, v(1));
app = diff(app, v(2));

pol_obj_2 = Polynomial(arb_basis);
pol2 = pol_obj_2.to_msspoly(); 

app = app + diff(pol2, v(2));

[vars, pows, coeffs] = decomp(app);

A_new = A_dec_var(pows, coeffs, vars, pol_obj, out_basis);

normm = norm(A - A_new)






