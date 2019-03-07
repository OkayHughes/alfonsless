%% problem description
%
% TODO: Write problem 
function out = find_FRS_alfonso(prob)

degree = prob.degree;

t = prob.t;
x = prob.x;
k = prob.k;


f = prob.f;
deg_f = ceil(max(arrayfun(@msspoly_degree, f))/2) * 2;

if isfield(prob, 'verbose')
  verbose = prob.verbose;
else
  verbose = false;
end

error_dynamics = isfield(prob, 'g');
if error_dynamics
    g = prob.g;
    num_gs = size(g, 2);
    deg_gs =  ceil(max(arrayfun(@msspoly_degree, g))/2) * 2;
end

X_bounds = prob.X_bounds;
X0_bounds = prob.X0_bounds;
K_bounds = prob.K_bounds;

T_min = 0;
T_max = prob.T;

%we make the following definitions:
% X = X_s
% Y = X x K
% Y0 = X0 x K
% Z = [0, T] x X x K

Y_vars = [x; k];
Y_bounds = [X_bounds;
            K_bounds];

Y0_bounds = [X0_bounds;
            K_bounds];

Z_vars = [t;x;k];
Z_bounds = [[T_min, T_max];
            X_bounds;
            K_bounds];


prog = AlfonsoSOSProgFekete;

prog.with_indeterminate([t;x;k]);

if verbose
  'Generating interpolant bases for Y, Y0'
end



time_index = 1;
first_sp_ind = 2;

[v, vcoeff] = prog.new_free_poly(Z_vars, degree);
[w, wcoeff, mon_w] = prog.new_free_poly(Y_vars, degree);

if error_dynamics
  q = msspoly(zeros([num_gs,1])) ;
  for qidx = 1:length(q)
      [q(qidx), ~] = prog.new_free_poly(Z_vars, degree) ;
  end
end

Lvf = diff(v, t) + diff(v, x) * f;
if error_dynamics
  Lvg = diff(v, x) * g;
end


%constraint 1
if verbose
  'Defining constraint 1'
end

if error_dynamics
    const1 = -(Lvf + ones(size(q))'*q);
    prog.sos_on_K(const1, Z_vars, Z_bounds,  degree + deg_f);
else
    const1 = -Lvf;
    prog.sos_on_K(const1, Z_vars, Z_bounds,  degree + deg_f);
end



if error_dynamics
  if verbose
  'Defining q constraints (2-4 in the paper)'
  end

  for g_ind=1:num_gs
      if verbose
        sprintf('Defining constraints for g_%d', g_ind)
      end
  %constraint 2
      const2 = Lvg(g_ind) + q(g_ind);
      prog.sos_on_K(const2, Z_vars, Z_bounds,  degree + deg_gs(g_ind));

  %constraint 3

      const3 = -Lvg(g_ind) + q(g_ind);
      prog.sos_on_K(const3, Z_vars, Z_bounds,  degree + deg_gs(g_ind));

  %constraint 4

      prog.sos_on_K(q(g_ind), Z_vars, Z_bounds, degree);
  end
end

if verbose
  'Defining constraint 5'
end

%constraint 5


const5 = -msubs(v, t, T_min);
prog.sos_on_K(const5, Y_vars, Y0_bounds, degree);


if verbose
  'Defining constraint 6'
end

%constraint 6


prog.sos_on_K(w, Y_vars, Y_bounds, degree);


if verbose
  'Defining constraint 7'
end

%constraint 7

const7 = w + v - 1;
prog.sos_on_K(const7, Z_vars, Z_bounds, degree);

if verbose
  'Defining cost function'
end

int_Y0 = boxMoments(Y_vars,Y_bounds(:,1),Y_bounds(:,2)) ;
obj = int_Y0(mon_w)'*(wcoeff);

if verbose
  'Running alfonso'
end

prog.problem_chars(2);
% run alfonso
exit
res = prog.minimize(obj);

res.polys;
out.v = res.polys(1);
out.w = res.polys(2);
out.qs = res.polys(3:end);
out.A = res.A;
out.b = res.b;
out.c = res.c;

if verbose
  'Done'
end

end

