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
    deg_g = ceil(max(arrayfun(@msspoly_degree, g))/2) * 2;
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
    int_params_Zg = AffineFeketeCube(Z_vars, degree + deg_g, Z_bounds);
end
if error_dynamics
    int_params_arr = [int_params_Zf, int_params_Zg, int_params_Zg, ...
                      int_params_Z, int_params_Y0, int_params_Y, int_params_Z];
else
    int_params_arr = [int_params_Zf, int_params_Y0, ...
                      int_params_Y, int_params_Z];
end

mon_basis_Y = int_params_Y.mon_basis;
mon_basis_Y0 = int_params_Y0.mon_basis;
mon_basis_Z = int_params_Z.mon_basis;
mon_basis_Zf = int_params_Zf.mon_basis;

if error_dynamics
  mon_basis_Zg = int_params_Zg.mon_basis;
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
    Q1 = int_params_Zf.mon_to_P0 * -Z_to_Zf * int_params_Z.P0_to_mon;
    prob_data.A = [V1, W1, Q1];
else
    prob_data.A = [V1, W1];
end

c1 = zeros(size(mon_basis_Zf.monomials));
prob_data.c = c1;

%constraint 2 (optional)

if error_dynamics
    L_g = liouville_operator_g(g, first_sp_ind, mon_basis_Z, mon_basis_Zg);
    V2 = int_params_Zg.mon_to_P0 * L_g * int_params_Z.P0_to_mon;
    W2 = zero_operator(mon_basis_Y, mon_basis_Zg);

    Z_to_Zg = monomial_to_monomial(mon_basis_Z, mon_basis_Zg);
    Q2 = int_params_Zg.mon_to_P0 * Z_to_Zg * int_params_Z.P0_to_mon;
    prob_data.A = [prob_data.A; 
                   [V2, W2, Q2]];
    c2 = zeros(size(mon_basis_Zg.monomials));
    prob_data.c = [prob_data.c;
                   c2];
end

%constraint 3 (optional)

if error_dynamics
    V3 = int_params_Zg.mon_to_P0 * -L_g * int_params_Z.P0_to_mon;
    W3 = zero_operator(mon_basis_Y, mon_basis_Zg);
    Z_to_Zg = monomial_to_monomial(mon_basis_Z, mon_basis_Zg);
    Q3 = int_params_Zg.mon_to_P0 * Z_to_Zg * int_params_Z.P0_to_mon;
    prob_data.A = [prob_data.A; 
                   [V3, W3, Q3]];
    c3 = zeros(size(mon_basis_Zg.monomials));
    prob_data.c = [prob_data.c;
                   c3];
end

%constraint 4 (optional)

if error_dynamics
    V4 = zero_operator(mon_basis_Z, mon_basis_Z);
    W4 = zero_operator(mon_basis_Y, mon_basis_Z);
    Q4 = eye(size(mon_basis_Z.monomials, 1));
    c4 = zeros(size(mon_basis_Z.monomials));
    prob_data.A = [prob_data.A; 
                   [V4, W4, Q4]];
    prob_data.c = [prob_data.c;
                   c4];
end

%constraint 5

apply_init_time = vector_partial_application(time_index, T_min, mon_basis_Z);
Z_to_Y0 = monomial_to_monomial(mon_basis_Z, mon_basis_Y0);
V5 = int_params_Y0.mon_to_P0 * Z_to_Y0 * -apply_init_time * int_params_Z.P0_to_mon;
W5 = zero_operator(mon_basis_Y, mon_basis_Y0);

if error_dynamics
    Q5 = zero_operator(mon_basis_Z, mon_basis_Y0);
    prob_data.A = [prob_data.A;
                   [V5, W5, Q5]];
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
    Q6 = zero_operator(mon_basis_Z, mon_basis_Y);
    prob_data.A = [prob_data.A;
                   [V6, W6, Q6]];
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
    Q7 = zero_operator(mon_basis_Z, mon_basis_Z);
    prob_data.A = [prob_data.A;
                   [V7, W7, Q7]];
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
    b3 = zeros(size(mon_basis_Z.monomials));
    prob_data.b = [prob_data.b;
                   b3];
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
    alfonso_q_vec = int_params_Z.P0_to_mon * results.y(int_params_Z.U+int_params_Y.U+1:end, 1);
    q_alfonso = alfonso_q_vec' * int_params_Z.mon_basis.monomials;

    out.q = q_alfonso;
end

out.w = w_alfonso;
out.v = v_alfonso;

end

