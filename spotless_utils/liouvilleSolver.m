% Liouville Solver with no input

function [w,v] = liouvilleSolver(t,x,f,hX,hXT,dl,T,degree)
%
%  t        -- 1-by-1 msspoly.
%  x{i}     -- n_i-by-1 free msspoly
%  f{i}     -- n_i-by-1 msspoly in x{i}
%  hX{i}    -- m_i-by-1 msspoly in x{i}
%  hXT{i}   -- k_i-by-1 msspoly in x{i}
%  dl{i}    -- function mapping n_i-by-1 msspoly's into n-by-1 double.
%  degree   -- positive scalar integer
%
%  Computes an outerapproximation to the BRS for:
%
%  xdot{i} = f_i(t,x{i}) + g_i(t,x{i})u,
%
%  on the time horizon [0,1] and set X, moving trajectories into
%  XT at t = 1, where:
%
%  X{i}  = { x | each hX{i}(x) >= 0 and sX{i} >= 0  }
%  XT{i} = { x | each hXT{i}(x) >= 0 }
%
%  dlambda(p) = int_X p(x) dx component-wise (lebesgue moments)
%
%  degree dictates the order of certain approximations.
%
%  w{i}(x) >= 1 is an outer approximation of x{i}(0) that could reach XT.
%          /some/ control law.

%% PROGRAM DEFINITION
prog = spotsosprog;

%% VARIABLES AND CONSTRAINTS
% Define the time constraint
hT  = t*(T-t);

% Make t and x indeterminate
prog = prog.withIndeterminate(t);
prog = prog.withIndeterminate(x) ;

% Create v polynomial
vmonom = monomials([t;x], 0:degree) ;
[prog, v, ~] = prog.newFreePoly(vmonom) ;

% Create w polynomial
wmonom = monomials(x, 0:degree) ;
[prog, w, wcoeff] = prog.newFreePoly(wmonom) ;

% Create subbed-in variables
v0 = subs(v, t, 0) ;
vT = subs(v, t, T) ;

% Liouville operator
Lv = diff(v,t) + diff(v,x)*f ;

% Create constraints 
% -Lv is SOS on[0,T] and X
prog = sosOnK(prog, -(Lv), [t;x], [hT; hX], degree) ;

% w - v0 - p - 1 is SOS on X
prog = sosOnK(prog, w - v0 - 1, x, hX, degree) ;

% vT - p is SOS on X_T
prog = sosOnK(prog, vT, x, hXT, degree) ;

% w is SOS on X
prog = sosOnK(prog, w, x, hX, degree) ;

% Create cost function
obj = dl(wmonom)' * wcoeff ;

%% SOLVER
% Set options
options = spot_sdp_default_options();
options.verbose = 1;

% Run program
sol = prog.minimize(obj, @spot_mosek, options);

% construct w{i}
w = sol.eval(w) ;
v = sol.eval(v) ;

end








