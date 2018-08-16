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

int_params_Y = AffineFeketeCube(Y_vars, degree, Y_bounds);
int_params_Y0 = AffineFeketeCube(Y_vars, degree, Y0_bounds);

if verbose
  'Generating interpolant bases for Z, Zf'
end

int_params_Z = AffineFeketeCube(Z_vars, degree, Z_bounds);
%Zf basis spans f_i*v
int_params_Zf = AffineFeketeCube(Z_vars, degree + deg_f, Z_bounds);
%Zf basis spans g_i*v
if error_dynamics
  if verbose
  'Generating interpolant bases for Zgs'
  end
  int_params_Zgs = repmat(int_params_Zf, 1, num_gs);
  for g_ind=1:num_gs
    deg_g_i = ceil(max(arrayfun(@msspoly_degree, g(:, g_ind)))/2) * 2;
    int_params_Zgs(g_ind) = AffineFeketeCube(Z_vars, degree + deg_g_i, Z_bounds);
  end
end

end