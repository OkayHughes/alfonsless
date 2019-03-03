
%% ALFONSO PROBLEM
% create interpolation points
function sol_alfonsless = SOS_fix_deg_exten(degree, mats)
dimension_param = 5;
variables = msspoly('x', dimension_param);

S_bounds = repmat([-1, 1], 5, 1);

%new alfonsless


dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));

prog = AlfonsoSOSProgFekete;
prog.with_indeterminate(variables);

[w, wcoeff, w_monomial] = prog.new_free_poly(variables, degree) ;
[v, vcoeff, v_monomial] = prog.new_free_poly(variables, degree) ;
[q, qcoeff, q_monomial] = prog.new_free_poly(variables, degree) ;
[r, rcoeff, r_monomial] = prog.new_free_poly(variables, degree) ;

for i=1:size(mats, 1)
	fprintf("defining const %d\n", i);
	pol = (mats{i, 1} * (wcoeff .* w_monomial) + mats{i, 2} * (vcoeff .* v_monomial) + mats{i, 3} * (qcoeff .* q_monomial) +  mats{i, 4} * (rcoeff .* r_monomial))';
	fprintf("adding to program\n")
	prog.sos_on_K(pol, variables, S_bounds, degree);
end

obj = dl(w_monomial)' * wcoeff;

fprintf("Running alfonsless\n");

sol_alfonsless = prog.minimize(obj);

end
