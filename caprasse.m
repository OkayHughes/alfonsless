%% problem description
% Solve the alfonso example problem "caprasse" using spotless. This does
% solve correctly.
%
%% user parameters
degree = 12 ;

%% automated from here
% variable t
t = msspoly('t',4) ;
t1 = t(1) ;
t2 = t(2) ;
t3 = t(3) ;
t4 = t(4) ;

% bounds on t
T_range = repmat([-0.5,0.5],length(t),1) ;
hT = (t - T_range(:,1)).*(T_range(:,2) - t) ;

% caprasse polynomial
f = -t1*t3^3 + 4*t2*t4*t3^2 + 4*t1*t3*t4^2 + 2*t2*t4^3 + 4*t1*t3 + ...
    +4*t3^2 - 10*t2*t4 - 10*t4^2 + 2 ;

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(t) ;

% create y variable
[prog,y] = prog.newFree(1) ;

% create SOS constraint
prog = sosOnK(prog, f - y, t, hT, degree) ;

% create cost function
obj = -y ;

%% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;

tic
sol = prog.minimize(obj, @spot_mosek, options);
toc