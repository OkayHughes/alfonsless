%% problem description
%
% TODO: Write problem 
function out = FRS_minimal_case(degree, name)

diary(name);
fprintf("min int_{T times X} v d mu s.t. v geq 0 on Z \n")
setup()
verbose = 1;
% options

% space parameters
T = 1 ;

%% set up advection vars
t = msspoly('t', 1);
x = msspoly('z', 2);
k = msspoly('k', 2) ;

T_min = 0;
T_max = T;
T_smol=0.1;
X_bounds = [-1, 1; -1, 1];
K_bounds = [-1, 1; -1, 1];


%we make the following definitions:
% X = X_s
% Y = X x K
% Y0 = X0 x K
% Z = [0, T] x X x K

Y_vars = [x; k];
Y_bounds = [X_bounds;
            K_bounds];

Z_vars = [t;x;k];
Z_bounds = [[T_min, T_max];
            X_bounds;
            K_bounds];

Z0_bounds = [[T_min, T_smol];
            X_bounds;
            K_bounds];


prog = AlfonsoSOSProgFekete;

prog.with_indeterminate([t;x;k]);


[v, vcoeff, v_mon] = prog.new_free_poly(Z_vars, degree);


prog.sos_on_K(v, Z_vars, Z_bounds, degree);

int_Z = boxMoments(Z_vars, Z_bounds(:,1), Z_bounds(:,2)) ;
obj = int_Z(v_mon)'*(vcoeff);

if verbose
  fprintf('Running alfonso\n')
end

prog.problem_chars(2);
res = prog.minimize(obj, struct());

res.polys;
out.v = res.polys(1);
out.A = res.A;
out.b = res.b;
out.c = res.c;

if verbose
  'Done'
end

diary off;
end

