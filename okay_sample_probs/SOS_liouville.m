%% problem description
%
% TODO: Write problem 
function SOS_liouville()

degree = 8 ;

t = msspoly('t');
x = msspoly('x');
variables = [t;x];
'dynamics:'
f = variables(2)^2
T_min = 0
T_max = 1

X_bounds = repmat([-1, 1], size(x, 1), 1)

intParams_w = AffineFeketeCube(x, degree, X_bounds);

X_T_bounds = [0.8, 1]

intParams_w_X_T = AffineFeketeCube(x, degree, X_T_bounds);

S_bounds = [[T_min, T_max];
            X_bounds];

intParams_v_f = AffineFeketeCube(variables, degree + ceil(msspoly_degree(f)/2) * 2, S_bounds);

intParams_arr = [intParams_v_f, intParams_w_X_T, intParams_w, intParams_w];

intParams_v = AffineFeketeCube(variables, degree, S_bounds);

gH_Params = gen_grad_params(intParams_arr);

v_to_v_f = monomial_to_monomial(intParams_v.mon_basis, intParams_v_f.mon_basis);
time_der_v = vector_derivative(1, intParams_v.mon_basis);
mult_mat = vector_poly_multiply(f, intParams_v.mon_basis, intParams_v_f.mon_basis);
space_der_v = vector_derivative(2, intParams_v.mon_basis);

A1 = -intParams_v_f.mon_to_P0 * liouville_operator(f, 1, 2, 2, intParams_v.mon_basis, intParams_v_f.mon_basis) * intParams_v.P0_to_mon;
A2 = zeros(intParams_v_f.U, intParams_w.U);

v_to_w_X_T = monomial_to_monomial(intParams_v.mon_basis, intParams_w_X_T.mon_basis);
apply_final_time = vector_partial_application(1, T_max, intParams_v.mon_basis);

B1 = intParams_w_X_T.mon_to_P0 * v_to_w_X_T * apply_final_time * intParams_v.P0_to_mon;
B2 = zeros(intParams_w.U);

v_to_w = monomial_to_monomial(intParams_v.mon_basis, intParams_w.mon_basis);
apply_init_time = vector_partial_application(1, T_min, intParams_v.mon_basis);

C1 = -intParams_w.mon_to_P0 * v_to_w * apply_init_time * intParams_v.P0_to_mon;
C2 = eye(intParams_w.U);

D1 = zeros(intParams_w.U, intParams_v.U);
D2 = eye(intParams_w.U);



const1 = zeros(intParams_v_f.U, 1);
const2 = zeros(intParams_w_X_T.U, 1);
const3 = ones(intParams_w.U, 1);
const4 = zeros(intParams_w.U, 1);

int1 = zeros(intParams_v.U, 1);
int2 = intParams_w.w;

probData.A = -[A1, A2;
     B1, B2;
     C1, C2;
     D1, D2]';
probData.c = -[const1;
      const2;
      const3;
      const4];
probData.b = -[int1;
      int2];


% make initial primal iterate
x0 = ones(sum(gH_Params.U_arr), 1);

[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),sum(gH_Params.U_arr),1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

hX = -(x-X_bounds(:, 1)).*(x-X_bounds(:, 2));
hXT = -(x - X_T_bounds(1)).*(x-X_T_bounds(2));
dl=boxMoments(x, X_bounds(:, 1), X_bounds(:, 2));

[w_spotless, v_spotless] = liouvilleSolver(t, x, f, hX, hXT, dl, T_max, degree);

%% PLOTTING POLYNOMIAL (RECOVERED) OUTPUT

spot_w_vec = msspoly_to_vector(w_spotless, intParams_w.mon_basis);
alfonso_w_vec = intParams_w.P0_to_mon * results.y(intParams_v.U+1:end, 1);
alfonso_v_vec = intParams_v.P0_to_mon * results.y(1:intParams_v.U, 1);
w_alfonso = alfonso_w_vec' * intParams_w.mon_basis.monomials;
v_alfonso = alfonso_v_vec' * intParams_v.mon_basis.monomials;
norm_v_subs = norm(dmsubs(v_alfonso, variables, intParams_v.pts')' - results.y(1:intParams_v.U, 1))
norm_v_f_subs = norm(intParams_v_f.mon_to_P0 * msspoly_to_vector(v_alfonso * f, intParams_v_f.mon_basis) - dmsubs(v_alfonso * f, variables, intParams_v_f.pts')')

%PLOTTING ALL CONSTRAINTS

v_to_v_f = monomial_to_monomial(intParams_v.mon_basis, intParams_v_f.mon_basis);
time_der_v = vector_derivative(1, intParams_v.mon_basis);
mult_mat = vector_poly_multiply(f, intParams_v.mon_basis, intParams_v_f.mon_basis);
space_der_v = vector_derivative(2, intParams_v.mon_basis);


constraint_1 = (intParams_v_f.P0_to_mon * A1 * results.y(1:intParams_v.U, 1))' * intParams_v_f.mon_basis.monomials;
constraint_1_alt = (-(v_to_v_f * time_der_v + mult_mat * space_der_v) * intParams_v.P0_to_mon * results.y(1:intParams_v.U, 1))' * intParams_v_f.mon_basis.monomials;

constraint_2 = (intParams_w_X_T.P0_to_mon * B1 *  results.y(1:intParams_v.U, 1))' * intParams_w.mon_basis.monomials;

v_to_w = monomial_to_monomial(intParams_v.mon_basis, intParams_w.mon_basis);
apply_init_time = vector_partial_application(1, T_min, intParams_v.mon_basis);

constraint_3 = (intParams_w.P0_to_mon * C1 * results.y(1:intParams_v.U, 1) + alfonso_w_vec)' * intParams_w.mon_basis.monomials + msspoly(-1);

constraint_4 = w_alfonso;

close all


figure('Name','v(T, x)') ; cla ; hold on ;
xlim(X_T_bounds);
ylim([-1, 1]);
pts = linspace(X_T_bounds(1, 1), X_T_bounds(1, 2), 100);
yvals = dmsubs(constraint_2, x, pts);
plot(pts, yvals)
plot(pts, zeros(size(pts)))

figure('Name', '-v(0, x) + w(x) - 1') ; cla ; hold on ;
xlim(X_bounds);
ylim([-1, 1]);
pts = linspace(X_bounds(1, 1), X_bounds(1, 2), 100);
yvals = dmsubs(constraint_3, x, pts);
plot(pts, yvals)
plot(pts, zeros(size(pts)))

figure('Name','w(x)') ; cla ; hold on ;
xlim(X_bounds);
ylim([-1, 1]);
pts = linspace(X_bounds(1, 1), X_bounds(1, 2), 100);
yvals = dmsubs(constraint_4, x, pts);
plot(pts, yvals)
plot(pts, zeros(size(pts)))

int_spotless = def_int_on_box(w_spotless, X_bounds(:, 1), X_bounds(:, 2), variables)
int_alfonso = def_int_on_box(w_alfonso, X_bounds(:, 1), X_bounds(:, 2), variables)

figure('Name','-(\mathcal{L}_f v)(x)') ; cla ; hold on ;
[grid_x, grid_y] = meshgrid(linspace(S_bounds(1, 1), S_bounds(1, 2), 50), linspace(S_bounds(2, 1), S_bounds(2, 2), 50));
flat = [grid_x(:), grid_y(:)];
const_1_vals = reshape(dmsubs(constraint_1, intParams_v_f.mon_basis.variables, flat')', 50, 50);
const_1_alt_vals = reshape(dmsubs(constraint_1_alt, intParams_v_f.mon_basis.variables, flat')', 50, 50);

COS(:,:,1) = zeros(50); % red
COS(:,:,2) = zeros(50); % green
COS(:,:,3) = zeros(50); % blue



surf_a = surf(grid_x, grid_y, const_1_vals, COS, 'FaceAlpha',0.5);
surf_c = surf(grid_x, grid_y, const_1_alt_vals, 'FaceAlpha',0);
surf_c.EdgeColor="red";
surf_b = surf(grid_x, grid_y, zeros(size(grid_x)), 'FaceAlpha',0);
surf_b.EdgeColor = 'white';




%[grid_x, grid_y] = meshgrid(linspace(-1, 1, 50), linspace(-1, 1, 50));
%flat = [grid_x(:), grid_y(:)];
x_grid = linspace(X_bounds(1, 1), X_bounds(1, 2), 50)';
alf_vals = dmsubs(w_alfonso, intParams_w.mon_basis.variables, x_grid')';
spot_vals = dmsubs(w_spotless, intParams_w.mon_basis.variables, x_grid')';



figure('Name','w_s(x) and w_a(x)') ; cla ; hold on ;

plot(x_grid, spot_vals,'LineWidth',1.5)
plot(x_grid,alf_vals,'--', 'LineWidth',1);

end
%surf(grid_x, grid_y, alf_vals, COA);
%surf(grid_x, grid_y, spot_vals, COS);
