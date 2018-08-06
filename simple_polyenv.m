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
fout = coeffs'*fmon ;

%% ALFONSO PROBLEM
% create interpolation points
intParams = ChebInterval(degree);
n   = intParams.n;
d   = intParams.d;
U   = intParams.U;
L   = intParams.L;
P   = intParams.P;
pts = intParams.pts;

% dimension of the weight polynomial space (should be dimension of d)
LWts = repmat(nchoosek(n+d-1,n),n,1);

% parameter object for cone gradient/hessian function
gH_Params.n = n;
gH_Params.d = d;
gH_Params.U = U;
gH_Params.numPolys = 2;

gH_Params.L     = L;
gH_Params.LWts  = LWts;
nu              = 2*(L+sum(LWts));
gH_Params.bnu	= nu+1;
gH_Params.P = P;

% create polynomial hT (g in the alfonso paper) to define space T = [-1,1]
wtVals = 1-pts.^2;
PWts{1}         = diag(sqrt(wtVals(:)))*P(:,1:LWts(1));
[PWts{1}, ~]    = qr(PWts{1}, 0);     
gH_Params.PWts = PWts;


% create SOS constraints by evaluating the polynomials f1 and f2 at the
% interpolation points
f1eval = msubs(f1,t,pts') ;
f2eval = msubs(f2,t,pts') ;

% create A,b,c matrices to define conic problem
probData.A = sparse(repmat(eye(U),1,2)) ;
probData.b = intParams.w ;
probData.c = -[f1eval(:); f2eval(:)] ;

% make initial primal iterate
x0 = ones(2*U,1) ;
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),2*U,1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

%% PLOTTING RESULTS
figure(1) ; cla ; hold on ;

% constraints
tvec = linspace(-1,1,100) ;
f1vals = msubs(f1,t,tvec) ;
f2vals = msubs(f2,t,tvec) ;
plot(tvec,f1vals,':','LineWidth',1)
plot(tvec,f2vals,':','LineWidth',1)

% spotless output
fvals = msubs(fout,t,tvec) ;
plot(tvec,fvals,'LineWidth',1.5)

% alfonso output:
yvals = -results.y ; % these are the values of the polynomial f at the points "pts"
plot(pts,yvals,'--','LineWidth',1.5)

xlabel('t')
ylabel('polynomial values')
legend('f1','f2','spotless','alfonso')