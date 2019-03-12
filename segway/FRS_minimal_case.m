%% problem description
%
% TODO: Write problem 
function out = FRS_minimal_case(degree, name)

diary(name);
setup()
verbose = 1;
% options

% initial conditions
x0 = -0.75 ;
y0 = 0 ;

% vehicle parameters
rcar = 0.2 ;
vmax = 0.5 ; % vmax = 2 usually
wmax = 1 ;

% space parameters
T = 1 ;
D = 1.5 ;

%% set up advection vars
t = msspoly('t', 1);
x = msspoly('z', 2);
z=x;
k = msspoly('k', 2) ;

kw = k(1) ;
kv = k(2) ;

%% set up spaces
X_bounds = repmat([-1, 1], size(z)) ;
X0_bounds = 1/D * [x0 - rcar, x0 + rcar;
                  y0 - rcar, y0 + rcar ];

K_bounds = repmat([-1,1], size(k));


T_min = 0;
T_max = T;

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


[v, vcoeff] = prog.new_free_poly(Z_vars, degree);
[w, wcoeff, mon_w] = prog.new_free_poly(Y_vars, degree);



if verbose
  'Defining constraint 7'
end

prog.sos_on_K(w, Y_vars, Y_bounds, degree);

% %constraint 7
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
%run alfonso
res = prog.minimize(obj, struct());

res.polys;
out.v = res.polys(1);
out.w = res.polys(2);
out.A = res.A;
out.b = res.b;
out.c = res.c;

if verbose
  'Done'
end

diary off;
end

