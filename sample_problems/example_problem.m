function example_problem()

fprintf("minimize \\int_{[-1, 2]} w(x) dx\ns.t. (u^3 + 3)w_x(u) \\geq 0 for u \\in [-1, 2]\n     w(u) \\geq 0 for u \\in [-1, 2]\n\n")

variables = msspoly('x');
fprintf("`variables = msspoly('x');`\nreturns:\n")
display(variables)
fprintf("`f = variables(1)^3+3;`\n returns:\n")
f=variables(1)^3+3;
display(f)

fprintf("`variables` defines the indeterminants in which decision polynomials will be defined, however in this case f is not a decision variable.\n If we set `variables = msspoly('x', 2)`, we could define f(x1, x2) = x1^2 + x2^2 with `f = variables(1)^2 + variables(2)^2;`\n\n ")

fprintf("`S_bounds = [-1, 2];`\n")
fprintf("Currently the only semialgebraic sets we can enforce positivity over are rectangular domains, i.e. [a_1, b_1] \\times [a_2, b_2] \\times ... \\times [a_n, b_n]\n `S_bounds` defines a rectangular domain.\n\n")

S_bounds = [-1, 2];


fprintf("`prog = AlfonsoSOSProgFekete;`\n")
fprintf("This creates an empty Alfonso program. \n")
fprintf("AlfonsoSOSProgFekete is a subclass of AlfonsoSOSProg which uses Fekete Cubes as the interpolant basis used in the Alfonso solver.\n")
fprintf("This is why all SOS constraints must be defined over rectangular domains for this subclass.\n\n")
prog = AlfonsoSOSProgFekete;

fprintf("`prog.with_indeterminate(variables);`\n")
fprintf("This line adds the indeterminates in `variables` to the program.\n\n")

prog.with_indeterminate(variables);

fprintf("`[w, wcoeff, wmonom] = prog.new_free_poly(variables, degree);`\n")
fprintf("`w` is a decision polynomial in variables `variables` (in this case in `x1`), with maximum degree 2. \n\n")

degree = 2;

[w, wcoeff, wmonom] = prog.new_free_poly(variables, degree) ;
fprintf("Note that `w` looks like:")
display(w)

fprintf("Note that this looks different from `f`, which we defined earlier:")
display(f)

fprintf("The msspoly vector `wcoeff` contains placeholder variables:")
wcoeff
fprintf("which correspond to the coefficients of the vector `wmonom`:");
wmonom
fprintf("Note that `wcoeff' * wmonom` is just w:");
wcoeff' * wmonom

fprintf("`bad_operator = diff(w, variables) * f;`\n")
fprintf("This calculates w_x1 using msspoly's differentiation function, and then multiplies by the polynomial `f`\n")

bad_operator = diff(w, variables) * f;

fprintf("Note that the actual variables over which we are optimizing correspond to entries of `wcoeff`,\n")
fprintf("and `f` contains none of the variables in `wcoeff`, so multiplying by f corresponds to a linear operation.\n\n")

fprintf("`prog.sos_on_K(poly, variables, S_bounds, degree);`\n")
fprintf("This creates a SOS constraint, ensuring that the polynomial `poly` in indeterminates `variables` is WSOS on the semialgebraic set defined by S_bounds.\n")
fprintf("The maximum degree of the SOS cone is given by the `degree` argument. In this case, we have multiplied by `f`, which is of degree 3, so we must include this in degree.\n\n ")
prog.sos_on_K(w, variables, S_bounds, degree);
prog.sos_on_K(bad_operator-1, variables, S_bounds, degree+ceil(msspoly_degree(f)/2) * 2);


fprintf("`dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));`\n")
fprintf("dl is a function which will take monomials with indeterminates in `variables` and return a real vector `v` with which `v'*wcoeff` calculates \\int_{S} w(x) dx\n;")
fprintf("In particular, v'*coeff is our")
dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));
fprintf("In particular, `v'*coeff` is the objective function we seek to optimize in this problem.")
obj = dl(wmonom)' * wcoeff ;

fprintf("`prog.minimize(obj)` formulates the problem and computes the Fekete Cubes that the Alfonsless solver uses to represent SOS cones\n")
fprintf("This will take by far the longest amount of time. Fekete Cubes will be saved to disk and reused later if the dimension and maximum degree match.\n")
fprintf("The function call returns a struct with as many decision polynomials as your problem calls for, as well as the real valued matrices that correspond to eqn (4) in the Alfonsless paper.\n")
opts.verbose = 0;
sol_alfonsless = prog.minimize(obj, opts);

fprintf("\n\nThe function call returns a struct with as many decision polynomials as your problem calls for,\n")
fprintf("as well as the real valued matrices that correspond to eqn (4) in the Alfonsless paper.\n")

fprintf("w_out = ")
sol_alfonsless.polys(1)
end
