
%% ALFONSO PROBLEM
% create interpolation points
function sol_alfonsless = OS_fix_deg_test(degree, A)
dimension_param = 5;
variables = msspoly('x', dimension_param);

S_bounds = repmat([-1, 1], 5, 1);



%new alfonsless


dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));

prog = AlfonsoSOSProgFekete;
prog.with_indeterminate(variables);

[w, wcoeff, w_monomial] = prog.new_free_poly(variables, degree) ;
pol = (A * (wcoeff .* w_monomial))' * ones(nchoosek(dimension_param + degree, dimension_param), 1)

prog.sos_on_K(pol, variables, S_bounds, degree);
prog.sos_on_K(w, variables, S_bounds, degree);
obj = dl(w_monomial)' * wcoeff ;

sol_alfonsless = prog.minimize(obj);

end
