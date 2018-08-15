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

int_params_Y = AffineFeketeCube(Y_vars, degree, Y_bounds);
int_params_Y0 = AffineFeketeCube(Y_vars, degree, Y0_bounds);

int_params_Z = AffineFeketeCube(Z_vars, degree, Z_bounds);
%Zf basis spans f_i*v
int_params_Zf = AffineFeketeCube(Z_vars, degree + deg_f, Z_bounds);
%Zf basis spans g_i*v
if error_dynamics
  int_params_Zgs = repmat(int_params_Zf, 1, num_gs);
  for g_ind=1:num_gs
    deg_g_i = ceil(max(arrayfun(@msspoly_degree, g(:, g_ind)))/2) * 2;
    int_params_Zgs(g_ind) = AffineFeketeCube(Z_vars, degree + deg_g_i, Z_bounds);
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

mon_basis_Y = int_params_Y.mon_basis;
mon_basis_Y0 = int_params_Y0.mon_basis;
mon_basis_Z = int_params_Z.mon_basis;
mon_basis_Zf = int_params_Zf.mon_basis;

if error_dynamics
  mon_basis_Zgs = repmat(int_params_Zgs(1).mon_basis, 1, num_gs);
  for g_ind=1:num_gs
    mon_basis_Zgs(g_ind) = int_params_Zgs(g_ind).mon_basis;
  end
end

%add intparams
g_h_params = gen_grad_params(int_params_arr);

time_index = 1;
first_sp_ind = 2;

L_f = liouville_operator(f, time_index, first_sp_ind, mon_basis_Z, mon_basis_Zf);
V1 = int_params_Zf.mon_to_P0 * -L_f * int_params_Z.P0_to_mon;
W1 = zero_operator(mon_basis_Y, mon_basis_Zf);


%constraint 1

if error_dynamics
    Z_to_Zf = monomial_to_monomial(mon_basis_Z, mon_basis_Zf);
    Q_arr_1 = repmat(int_params_Zf.mon_to_P0 * -Z_to_Zf * int_params_Z.P0_to_mon, 1, num_gs);
    prob_data.A = [V1, W1, Q_arr_1];
else
    prob_data.A = [V1, W1];
end

c1 = zeros(size(mon_basis_Zf.monomials));
prob_data.c = c1;



if error_dynamics

  for g_ind=1:num_gs
  %constraint 2
      int_params_Zg = int_params_Zgs(g_ind);
      mon_basis_Zg = mon_basis_Zgs(g_ind);
      L_g = liouville_operator_g(g(:, g_ind), first_sp_ind, mon_basis_Z, mon_basis_Zg);
      V2 = int_params_Zg.mon_to_P0 * L_g * int_params_Z.P0_to_mon;
      W2 = zero_operator(mon_basis_Y, mon_basis_Zg);

      Z_to_Zg = monomial_to_monomial(mon_basis_Z, mon_basis_Zg);
      Q_i_2 = int_params_Zg.mon_to_P0 * Z_to_Zg * int_params_Z.P0_to_mon;
      Q_zeros_2 = zero_operator(mon_basis_Z, mon_basis_Zg);

      Q_arr_2 = [repmat(Q_zeros_2, 1, g_ind-1), Q_i_2, repmat(Q_zeros_2, 1, num_gs - g_ind)];
      prob_data.A = [prob_data.A; 
                     [V2, W2, Q_arr_2]];
      c2 = zeros(size(mon_basis_Zg.monomials));
      prob_data.c = [prob_data.c;
                     c2];

  %constraint 3

      V3 = int_params_Zg.mon_to_P0 * -L_g * int_params_Z.P0_to_mon;
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

%constraint 5

apply_init_time = vector_partial_application(time_index, T_min, mon_basis_Z);
Z_to_Y0 = monomial_to_monomial(mon_basis_Z, mon_basis_Y0);
V5 = int_params_Y0.mon_to_P0 * Z_to_Y0 * -apply_init_time * int_params_Z.P0_to_mon;
W5 = zero_operator(mon_basis_Y, mon_basis_Y0);

if error_dynamics
    Q_zeros_5 = zero_operator(mon_basis_Z, mon_basis_Y0);
    Q_arr_5 = repmat(Q_zeros_5, 1, num_gs);
    prob_data.A = [prob_data.A;
                   [V5, W5, Q_arr_5]];
else
    prob_data.A = [prob_data.A;
                   [V5, W5]];
end

c5 = zeros(size(mon_basis_Y0.monomials));

prob_data.c = [prob_data.c;
               c5];

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

%constraint 7

Y_to_Z = monomial_to_monomial(mon_basis_Y, mon_basis_Z);

V7 = eye(size(mon_basis_Z.monomials, 1));
W7 = int_params_Z.mon_to_P0 * Y_to_Z * int_params_Y.P0_to_mon;

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


b1 = zeros(size(mon_basis_Z.monomials));
prob_data.b = b1;

b2 = -int_params_Y.w;
prob_data.b = [prob_data.b;
               b2];

if error_dynamics
  for g_ind=1:num_gs
    b_g_i = zeros(size(mon_basis_Z.monomials));
    prob_data.b = [prob_data.b;
                   b_g_i];
  end
end

% -Af + s = -c => Af >= c
prob_data.A = -transpose(prob_data.A);
prob_data.c = -prob_data.c;


% make initial primal iterate
x0 = ones(sum(g_h_params.U_arr), 1);
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, g_h_params);

rP = max((1+abs(prob_data.b))./(1+abs(prob_data.A*x0)));
rD = max((1+abs(g0))./(1+abs(prob_data.c)));
x0 = repmat(sqrt(rP*rD),sum(g_h_params.U_arr),1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(prob_data, x0, @alfonso_grad_and_hess, g_h_params, opts);

alfonso_w_vec = int_params_Y.P0_to_mon * results.y(int_params_Z.U+1:int_params_Z.U+int_params_Y.U, 1);
alfonso_v_vec = int_params_Z.P0_to_mon * results.y(1:int_params_Z.U, 1);
w_alfonso = alfonso_w_vec' * int_params_Y.mon_basis.monomials;
v_alfonso = alfonso_v_vec' * int_params_Z.mon_basis.monomials;


if error_dynamics
    out.qs = msspoly(zeros(num_gs, 1));
    for g_ind=1:num_gs
      low_ind = int_params_Z.U + int_params_Y.U + int_params_Z.U * (g_ind-1)+ 1;
      high_ind = int_params_Z.U + int_params_Y.U + int_params_Z.U * (g_ind);
      alfonso_q_vec = int_params_Z.P0_to_mon * results.y(low_ind:high_ind, 1);
      q_alfonso = alfonso_q_vec' * int_params_Z.mon_basis.monomials;

      out.qs(g_ind) = q_alfonso;
    end
end

out.w = w_alfonso;
out.v = v_alfonso;

end

