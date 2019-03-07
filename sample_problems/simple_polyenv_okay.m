%% problem description
% This program does actually work.
%
%   min  int f(t)dt
%    f    T
%
%   s.t. f >= 0.5 - t^2
%        f >= 0.5 - (t-0.25)^2
%
%% user parameters
degree = 30 ; % degree of f

%% SPOTLESS PROBLEM
% polynomial variable
t = msspoly('t',1) ;

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(t);

% bounds on t
T_range = [-1,1];
hT = (t - T_range(:,1)).*(T_range(:,2) - t) ;

% create constraint polynomials
f1 = 0.5 - t^2 ;
f2 = 0.5 - (t-0.25)^2 ;


% create decision variable
fmon = monomials(t,0:degree) ;
[prog,f,fcoeff] = prog.newFreePoly(fmon) ;

% create SOS constraints
prog = sosOnK(prog, f - f1, t, hT, degree) ;
prog = sosOnK(prog, f - f2, t, hT, degree) ;

% create cost function
int_T = boxMoments(t,T_range(:,1),T_range(:,2)) ;
obj = int_T(fmon)'*fcoeff ;

% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;

tic
sol = prog.minimize(obj, @spot_mosek, options);
toc

% recover f
coeffs = sol.eval(fcoeff) ;

fspotless = coeffs'*fmon ;

cost_spotless = dmsubs(obj, fcoeff, dmsubs(coeffs, fcoeff, ones(size(coeffs))))


prog = AlfonsoSOSProgFekete;
prog.with_indeterminate(t);


% create decision variable
[f,fcoeff] = prog.new_free_poly(t, degree) ;

% create SOS constraints
prog.sos_on_K(f - f1, T_range, degree) ;
prog.sos_on_K(f - f2, T_range, degree) ;

% create cost function
int_T = boxMoments(t,T_range(:,1),T_range(:,2)) ;
obj = int_T(fmon)'*fcoeff ;

res = prog.minimize(obj);
cost_alfonso = res.cost

falfonso = res.polys(1);

l2_dist = l2_dist_on_box(fspotless, falfonso, T_range(:, 1), T_range(:, 2), t)



%% PLOTTING RESULTS
figure(1) ; cla ; hold on ;

% % constraints
% tvec = linspace(-1,1,100) ;
% f1vals = msubs(f1,t,tvec) ;
% f2vals = msubs(f2,t,tvec) ;
% plot(tvec,f1vals,':','LineWidth',1)
% plot(tvec,f2vals,':','LineWidth',1)

% spotless output
fs = msubs(fspotless,t,tvec);
plot(tvec,fs,'LineWidth',1.5)

fa = msubs(falfonso,t,tvec);
plot(tvec,fa,'LineWidth',1.5)

xlabel('t')
ylabel('polynomial values')
legend('f1','f2','spotless','alfonso')