%% problem description
%
% TODO: Write problem description
degree = 8 ;

t = msspoly('t');
x = msspoly('x');
variables = [t;x];
'dynamics:'
f = variables(2)^2
T_min = 0
T_max = 1

X_bounds = repmat([-1, 1], size(x, 1), 1)

intParams_w = FeketeCube(size(x, 1), degree/2, x);
intParams_w = scale_fekete_cube(intParams_w, X_bounds);

X_T_bounds = [-0.1, 0.1]

intParams_w_X_T = FeketeCube(size(x, 1), degree/2, x);
intParams_w_X_T = scale_fekete_cube(intParams_w_X_T, X_T_bounds);

S_bounds = [[T_min, T_max];
            repmat([-1, 1], size(x, 1), 1)];

intParams_v_f = FeketeCube(size(variables, 1), degree/2 + ceil(msspoly_degree(f)/2), variables); 
intParams_v_f = scale_fekete_cube(intParams_v_f, S_bounds);
intParams_arr = [intParams_v_f, intParams_w_X_T, intParams_w, intParams_w];


intParams_v = FeketeCube(size(variables, 1), degree/2, variables);
intParams_v = scale_fekete_cube(intParams_v, S_bounds);


numPolys = 4;

gH_Params.numPolys = numPolys;
gH_Params.n_arr = zeros(numPolys, 1);
gH_Params.d_arr = zeros(numPolys, 1);
gH_Params.U_arr = zeros(numPolys, 1);
gH_Params.L_arr = zeros(numPolys, 1);
gH_Params.LWts_cell = cell(numPolys, 1);
gH_Params.P_cell = cell(numPolys, 1);
gH_Params.PWts_cell = cell(numPolys, 1);
gH_Params.bnu = 0;

for i=1:numPolys
    intParams_i = intParams_arr(i);
    n   = intParams_i.n;
    d   = intParams_i.d;
    U   = intParams_i.U;
    L   = intParams_i.L;
    P   = intParams_i.P0;
    pts = intParams_i.pts;
    bounds = intParams_i.bounds;

    % dimension of the weight polynomial space (should be dimension of d)
    LWts = repmat(nchoosek(n+d-1,n),n,1);
    % parameter object for cone gradient/hessian function
    gH_Params.n_arr(i)      = n;
    gH_Params.d_arr(i)      = d;
    gH_Params.U_arr(i)      = U;
    gH_Params.L_arr(i)      = L;
    gH_Params.LWts_cell{i}  = LWts;
    nu                      = L+sum(LWts) ;
    gH_Params.bnu           = gH_Params.bnu + nu;
    gH_Params.P_cell{i}     = P;

    lb = bounds(:, 1);
    ub = bounds(:, 2);
    % create polynomial hT (g in the alfonso paper) to define space T = [-1,1]^2

    wtVals  = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);

    PWts = cell(n,1);
    for j = 1:n
        PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
        [PWts{j}, ~]    = qr(PWts{j}, 0);
        % associated positive semidefinite cone constraints: 
        % PWts{j}'*diag(x_1)*PWts{j} >= 0,
        % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
    end
    gH_Params.PWts_cell{i} = PWts;
end
gH_Params.bnu = gH_Params.bnu + 1;

v_to_v_f = monomial_to_monomial(intParams_v.mon_basis, intParams_v_f.mon_basis);
time_der_v = vector_derivative(1, intParams_v.mon_basis);
mult_mat = vector_poly_multiply(f, intParams_v.mon_basis, intParams_v_f.mon_basis);
space_der_v = vector_derivative(2, intParams_v.mon_basis);

A1 = intParams_v_f.mon_to_P0 * (v_to_v_f * time_der_v + mult_mat * space_der_v) * intParams_v.P0_to_mon;
A2 = zeros(intParams_v_f.U, intParams_w.U);

v_to_w_X_T = monomial_to_monomial(intParams_v.mon_basis, intParams_w_X_T.mon_basis);
apply_final_time = vector_partial_application(1, T_max, intParams_v.mon_basis);

B1 = intParams_w_X_T.mon_to_P0 * v_to_w_X_T * apply_final_time * intParams_v.P0_to_mon;
B2 = zeros(intParams_w.U);

v_to_w = monomial_to_monomial(intParams_v.mon_basis, intParams_w.mon_basis);
apply_init_time = vector_partial_application(1, T_min, intParams_v.mon_basis);

C1 = intParams_w.mon_to_P0 * v_to_w * apply_init_time * intParams_v.P0_to_mon;
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
hXT = -(t - X_T_bounds(1))*(t-X_T_bounds(2));
dl=boxMoments(x, X_bounds(:, 1), X_bounds(:, 2));

[w_spotless, spotless_v] = liouvilleSolver(t, x, f, hX, hXT, dl, T_max, degree);

%% PLOTTING POLYNOMIAL (RECOVERED) OUTPUT

spot_poly_vec = msspoly_to_vector(w_spotless, intParams_w.mon_basis);
alfonso_poly_vec = intParams_w.P0_to_mon * results.y(intParams_v.U+1:end, 1);
falfonso = alfonso_poly_vec' * intParams_w.mon_basis.monomials;


%[grid_x, grid_y] = meshgrid(linspace(-1, 1, 50), linspace(-1, 1, 50));
%flat = [grid_x(:), grid_y(:)];
x_grid = linspace(0, 1, 50)';
alf_vals = dmsubs(falfonso, intParams_w.mon_basis.variables, x_grid')';
spot_vals = dmsubs(w_spotless, intParams_w.mon_basis.variables, x_grid')';


figure(1) ; cla ; hold on ;

plot(x_grid, spot_vals,'LineWidth',1.5)
size(x_grid)
size(alf_vals)
plot(x_grid,alf_vals,'LineWidth',1);
%surf(grid_x, grid_y, alf_vals, COA);
%surf(grid_x, grid_y, spot_vals, COS);
