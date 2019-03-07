%% problem description
% This problem does not always solve correctly. The point "p" must be one
% of the interpolation points for spotless and alfonso to yield the same
% solution. The stuff defining alfonso's A, b, and c matrices is broken.
%
%   min   int f(t)dt
%    f     T
%
%   s.t.  f >= 0
%         f(p) >= c
%
%% user parameters
degree = 12 ;

% point constraint in the problem description
p = 0.2; 
c = 1.5 ;

%% SPOTLESS PROBLEM
% polynomial variable
t = msspoly('t', 1);

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(t) ;

% bounds on t
T_range = [-1,1] ;
hT = (t - T_range(:,1)).*(T_range(:,2) - t) ;

% create decision variables
mon = monomials(t,0:degree) ;
[prog,f,fcoeff] = prog.newFreePoly(mon) ;

% create SOS constraints
fcon = msubs(f,t,p) ;
prog = sosOnK(prog, f, t, hT, degree) ;
prog = sosOnK(prog, fcon - c, t, hT, degree) ;

% create cost function
int_T = boxMoments(t,T_range(:,1),T_range(:,2)) ;
obj = int_T(mon)'*(fcoeff);

% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;
sol = prog.minimize(obj, @spot_mosek, options);

% recover f
fspotless = sol.eval(fcoeff)'*mon ;



prog = AlfonsoSOSProgFekete ;
prog.with_indeterminate(t) ;

% bounds on t
T_range = [-1,1] ;

% create decision variables
[f,fcoeff, mon] = prog.new_free_poly(t, degree) ;

% create SOS constraints
fcon = msubs(f,t,p) ;
prog.sos_on_K(f, T_range, degree) ;
prog.sos_on_K(fcon - c, T_range, degree) ;

% create cost function
int_T = boxMoments(t,T_range(:,1),T_range(:,2)) ;
obj = int_T(mon)'*(fcoeff);

% run solver
sol = prog.minimize(obj);

% recover f
falfonso = sol.polys(1);



%% PLOTTING RESULTS
figure(1) ; cla ; hold on ;
xlim([-1, 1])
ylim([-1, 2])

tvec = linspace(-1,1,500) ;

% spotless output
fvals = msubs(fspotless,t,tvec) ;
plot(tvec,fvals,'LineWidth',1.5)

% alfonso output:
%yvals = -results.y ; % these are the values of the polynomial f at the points "pts"

monvals = dmsubs(falfonso,t,tvec);
scatter(tvec,monvals,'LineWidth',1);
%results.y)
% figure(1) ; cla ; hold on ;

% xlim([-1, 1])
% ylim([-1, 2])

