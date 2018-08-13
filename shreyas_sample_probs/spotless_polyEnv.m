%% problem description
% This does work correctly!
%
% solve the polynomial envelope example problem from alfonso using spotless
%
%% user parameters
n = 1 ;
degree = 100 ;
numPolys = 2 ;
degPolys = 5 ;

seed = 2017 ;

%% automated from here
% seed the RNG
rng(seed, 'twister') ;

% variable t
t = msspoly('t',n) ;

% bounds on t
T_range = [-rand(n,1),rand(n,1)] ;
hT = (t - T_range(:,1)).*(T_range(:,2) - t) ;

% create constraint polynomials
LDegs = nchoosek(n+degPolys, n);
coeffs = randi([-9,9], LDegs, numPolys);
mmon = monomials(t,0:degPolys) ;
f = coeffs'*mmon ;

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(t) ;

% create y polynomial
ymon = monomials(t,0:degree) ;
[prog,y,ycoeff] = prog.newFreePoly(ymon) ;

% create SOS constraints
for idx = 1:size(f,1)
    fi = f(idx) ;
    prog = sosOnK(prog, y - f, t, hT, degree) ;
end

% create cost function
int_T = boxMoments(t,T_range(:,1),T_range(:,2)) ;
obj = int_T(ymon)'*ycoeff ;

%% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;

tic
sol = prog.minimize(obj, @spot_mosek, options);
toc