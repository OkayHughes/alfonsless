
%% ALFONSO PROBLEM
% create interpolation points
function sol_alfonsless = SOS_fix_deg_exten(num_polys, degree, mats)
dimension_param = 5;
variables = msspoly('x', dimension_param);

S_bounds = repmat([-1, 1], 5, 1);

%new alfonsless


dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));

prog = AlfonsoSOSProgFekete;
prog.with_indeterminate(variables);

polys = msspoly();

coeffs = msspoly();

mons = msspoly();

for i=1:num_polys
    [poly, coeff, mon] = prog.new_free_poly(variables, degree) ;
    polys = [polys; poly];
    coeffs = [coeffs; coeff'];
    mons = [mons; mon'];
end

for i=1:size(mats, 1)
    pol = msspoly(0);
	fprintf("defining const %d\n", i);
    for j=1:num_polys
        pol = pol + ((mats{i, j}' .* coeffs(j, :)')' * mons(j, :)');
    end
    fprintf("adding to program\n")
	prog.sos_on_K(pol, variables, S_bounds, degree);
end

obj = dl(mons(1, :)')' * coeffs(1, :)';

fprintf("Running alfonsless\n");

sol_alfonsless = prog.minimize(obj);

end
