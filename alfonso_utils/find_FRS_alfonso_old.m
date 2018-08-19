%% problem description
%
% TODO: Write problem 
function out = find_FRS_alfonso_old(prob)

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

if verbose
  'Generating interpolant bases for Y, Y0'
end

int_params_Y = AffineFeketeBasis(Y_vars, degree, Y_bounds);
int_params_Y0 = AffineFeketeBasis(Y_vars, degree, Y0_bounds);

if verbose
  'Generating interpolant bases for Z, Zf'
end

int_params_Z = AffineFeketeBasis(Z_vars, degree, Z_bounds);
%Zf basis spans f_i*v
int_params_Zf = AffineFeketeBasis(Z_vars, degree + deg_f, Z_bounds);
%Zf basis spans g_i*v
if error_dynamics
  if verbose
  'Generating interpolant bases for Zgs'
  end
  int_params_Zgs = repmat(int_params_Zf, 1, num_gs);
  for g_ind=1:num_gs
    deg_g_i = ceil(max(arrayfun(@msspoly_degree, g(:, g_ind)))/2) * 2;
    int_params_Zgs(g_ind) = AffineFeketeBasis(Z_vars, degree + deg_g_i, Z_bounds);
  end
end

if error_dynamics
    int_params_arr = int_params_Zf;
    for g_ind=1:num_gs
      int_params_arr = [int_params_arr, int_params_Zgs(g_ind), int_params_Zgs(g_ind), int_params_Z];
    end
    int_params_arr = [int_params_arr, int_params_Y0, int_params_Y, int_params_Z];
else
    int_params_arr = [int_params_Zf, int_params_Y0, ...
                      int_params_Y, int_params_Z];
end

mon_basis_Y = MonomialBasis(int_params_Y.variables, int_params_Y.d);
mon_basis_Z = MonomialBasis(int_params_Z.variables, int_params_Z.d);
mon_basis_Zf = MonomialBasis(int_params_Zf.variables, int_params_Zf.d);

if error_dynamics
  mon_basis_Zgs = repmat(MonomialBasis(int_params_Zgs(1).variables, int_params_Zgs(1).d), 1, num_gs);
  for g_ind=1:num_gs
    mon_basis_Zgs(g_ind) = MonomialBasis(int_params_Zgs(g_ind).variables, int_params_Zgs(g_ind).d);
  end
end

if verbose
  'Generating g_h_params'
end
%add intparams
g_h_params = gen_grad_params(int_params_arr);

time_index = 1;
first_sp_ind = 2;

[Zf_P0_to_mon, Zf_mon_to_P0] = int_params_Zf.inter_to_mon(mon_basis_Zf);
[Z_P0_to_mon, Z_mon_to_P0] = int_params_Z.inter_to_mon(mon_basis_Z);



L_f = liouville_operator(f, time_index, first_sp_ind, mon_basis_Z, mon_basis_Zf);
V1 = Zf_mon_to_P0 * -L_f * Z_P0_to_mon;
W1 = zero_operator(mon_basis_Y, mon_basis_Zf);


%constraint 1
if verbose
  'Defining constraint 1'
end

if error_dynamics
    Z_to_Zf = monomial_to_monomial(mon_basis_Z, mon_basis_Zf);
    Q_arr_1 = repmat(Zf_mon_to_P0 * -Z_to_Zf * Z_P0_to_mon, 1, num_gs);
    prob_data.A = [V1, W1, Q_arr_1];
else
    prob_data.A = [V1, W1];
end

c1 = zeros(size(mon_basis_Zf.monomials));
prob_data.c = c1;



if error_dynamics
  if verbose
  'Defining q constraints (2-4 in the paper)'
  end

  for g_ind=1:num_gs
      if verbose
       sprintf('Defining constraints for g_%d', g_ind)
      end
  %constraint 2
      int_params_Zg = int_params_Zgs(g_ind);
      mon_basis_Zg = mon_basis_Zgs(g_ind);
      [Zg_P0_to_mon, Zg_mon_to_P0] = int_params_Zg.inter_to_mon(mon_basis_Zg);
      
      L_g = liouville_operator_g(g(:, g_ind), first_sp_ind, mon_basis_Z, mon_basis_Zg);
      V2 = Zg_mon_to_P0 * L_g * Z_P0_to_mon;
      W2 = zero_operator(mon_basis_Y, mon_basis_Zg);

      Z_to_Zg = monomial_to_monomial(mon_basis_Z, mon_basis_Zg);
      Q_i_2 = Zg_mon_to_P0 * Z_to_Zg * Z_P0_to_mon;
      Q_zeros_2 = zero_operator(mon_basis_Z, mon_basis_Zg);

      Q_arr_2 = [repmat(Q_zeros_2, 1, g_ind-1), Q_i_2, repmat(Q_zeros_2, 1, num_gs - g_ind)];
      prob_data.A = [prob_data.A; 
                     [V2, W2, Q_arr_2]];
      c2 = zeros(size(mon_basis_Zg.monomials));
      prob_data.c = [prob_data.c;
                     c2];

  %constraint 3

      V3 = Zg_mon_to_P0 * -L_g * Z_P0_to_mon;
      W3 = zero_operator(mon_basis_Y, mon_basis_Zg);
      Q_i_3 = Q_i_2;
      Q_zeros_3 = Q_zeros_2;
      Q_arr_3 = [repmat(Q_zeros_3, 1, g_ind-1), Q_i_3, repmat(Q_zeros_3, 1, num_gs - g_ind)];
      prob_data.A = [prob_data.A; 
                     [V3, W3, Q_arr_3]];
      c3 = zeros(size(mon_basis_Zg.monomials));
      prob_data.c = [prob_data.c;
                     c3];

  %constraint 4

      V4 = zero_operator(mon_basis_Z, mon_basis_Z);
      W4 = zero_operator(mon_basis_Y, mon_basis_Z);
      Q_i_4 = eye(size(mon_basis_Z.monomials, 1));
      Q_zeros_4 = zero_operator(mon_basis_Z, mon_basis_Z);
      Q_arr_4 = [repmat(Q_zeros_4, 1, g_ind-1), Q_i_4, repmat(Q_zeros_4, 1, num_gs - g_ind)];

      c4 = zeros(size(mon_basis_Z.monomials));
      prob_data.A = [prob_data.A; 
                     [V4, W4, Q_arr_4]];
      prob_data.c = [prob_data.c;
                     c4];
  end
end

if verbose
  'Defining constraint 5'
end

%constraint 5

[Y0_P0_to_mon, Y0_mon_to_P0] = int_params_Y0.inter_to_mon(mon_basis_Y);

apply_init_time = vector_partial_application(time_index, T_min, mon_basis_Z);
Z_to_Y0 = monomial_to_monomial(mon_basis_Z, mon_basis_Y);
V5 = Y0_mon_to_P0 * Z_to_Y0 * -apply_init_time * Z_P0_to_mon;
W5 = zero_operator(mon_basis_Y, mon_basis_Y);

if error_dynamics
    Q_zeros_5 = zero_operator(mon_basis_Z, mon_basis_Y);
    Q_arr_5 = repmat(Q_zeros_5, 1, num_gs);
    prob_data.A = [prob_data.A;
                   [V5, W5, Q_arr_5]];
else
    prob_data.A = [prob_data.A;
                   [V5, W5]];
end

c5 = zeros(size(mon_basis_Y.monomials));

prob_data.c = [prob_data.c;
               c5];

if verbose
  'Defining constraint 6'
end

%constraint 6

V6 = zero_operator(mon_basis_Z, mon_basis_Y);
W6 = eye(size(mon_basis_Y.monomials, 1));

if error_dynamics
    Q_zeros_6 = zero_operator(mon_basis_Z, mon_basis_Y);
    Q_arr_6 = repmat(Q_zeros_6, 1, num_gs);
    prob_data.A = [prob_data.A;
                   [V6, W6, Q_arr_6]];
else
    prob_data.A = [prob_data.A;
                   [V6, W6]];
end

c6 = zeros(size(mon_basis_Y.monomials));

prob_data.c = [prob_data.c;
               c6];

if verbose
  'Defining constraint 7'
end

%constraint 7

[Y_P0_to_mon, Y_mon_to_P0] = int_params_Y.inter_to_mon(mon_basis_Y);

Y_to_Z = monomial_to_monomial(mon_basis_Y, mon_basis_Z);

V7 = eye(size(mon_basis_Z.monomials, 1));
W7 = Z_mon_to_P0 * Y_to_Z * Y_P0_to_mon;

if error_dynamics
    Q_zeros_7 = zero_operator(mon_basis_Z, mon_basis_Z);
    Q_arr_7 = repmat(Q_zeros_7, 1, num_gs);
    prob_data.A = [prob_data.A;
                   [V7, W7, Q_arr_7]];
else
    prob_data.A = [prob_data.A;
                   [V7, W7]];
end

c7 = ones(size(mon_basis_Z.monomials));

prob_data.c = [prob_data.c;
               c7];

if verbose
  'Defining cost function'
end

b1 = zeros(size(mon_basis_Z.monomials));
prob_data.b = b1;
size(prob_data.b)

b2 = -int_params_Y.w;
prob_data.b = [prob_data.b;
               b2];
size(prob_data.b)
if error_dynamics
  for g_ind=1:num_gs
    b_g_i = zeros(size(mon_basis_Z.monomials));
    prob_data.b = [prob_data.b;
                   b_g_i];
    size(prob_data.b)
  end
end
% -Af + s = -c => Af >= c
prob_data.A = -transpose(prob_data.A);
prob_data.c = -prob_data.c;

if verbose
  'Creating primal iterate'
end

% make initial primal iterate
x0 = ones(sum(g_h_params.U_arr), 1);
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, g_h_params);

rP = max((1+abs(prob_data.b))./(1+abs(prob_data.A*x0)));
rD = max((1+abs(g0))./(1+abs(prob_data.c)));
x0 = repmat(sqrt(rP*rD),sum(g_h_params.U_arr),1);

if verbose
  'Running alfonso'
end

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(prob_data, x0, @alfonso_grad_and_hess, g_h_params, opts);

if verbose
  'Extracting results'
end

alfonso_w_vec = Y_P0_to_mon * results.y(int_params_Z.U+1:int_params_Z.U+int_params_Y.U, 1);
alfonso_v_vec = Z_P0_to_mon * results.y(1:int_params_Z.U, 1);
w = alfonso_w_vec' * mon_basis_Y.monomials;
v = alfonso_v_vec' * mon_basis_Z.monomials;


if error_dynamics
    out.qs = msspoly(zeros(num_gs, 1));
    for g_ind=1:num_gs
      low_ind = int_params_Z.U + int_params_Y.U + int_params_Z.U * (g_ind-1)+ 1;
      high_ind = int_params_Z.U + int_params_Y.U + int_params_Z.U * (g_ind);
      alfonso_q_vec = Z_P0_to_mon * results.y(low_ind:high_ind, 1);
      q = alfonso_q_vec' * mon_basis_Z.monomials;

      out.qs(g_ind) = q;
    end
end

out.w = w;
out.v = v;

out.A = prob_data.A;
out.b = prob_data.b;
out.c = prob_data.c;

if verbose
  'Done'
end

end
