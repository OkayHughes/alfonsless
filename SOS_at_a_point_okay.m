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
intParams = FeketeCube(1,degree) ;

% point constraint in the problem description
p = 0.5; 
c = 1.5 ;

%% SPOTLESS PROBLEM
% polynomial variable
t = intParams.mon_basis.variables ;

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
fmosek = sol.eval(fcoeff)'*mon ;

%% ALFONSO PROBLEM
% create interpolation points
n   = intParams.n;
d   = intParams.d;
U   = intParams.U;
L   = intParams.L;
P   = intParams.P0;
pts = intParams.pts;

% dimension of the weight polynomial space (should be dimension of d)
LWts = repmat(nchoosek(n+d-1,n),n,1);

% parameter object for cone gradient/hessian function
gH_Params.n = n;
gH_Params.d = d;
gH_Params.U = U;
gH_Params.numPolys = 2; % there are two SOS constraints

gH_Params.L     = L;
gH_Params.LWts  = LWts;
nu              = gH_Params.numPolys*(L+sum(LWts)) ;
gH_Params.bnu   = nu+1 ;

gH_Params.P = P;

% create polynomial hT (g in the alfonso paper) to define space T = [-1,1]
wtVals = 1-pts.^2;
PWts{1}         = diag(sqrt(wtVals(:)))*P(:,1:LWts(1));
[PWts{1}, ~]    = qr(PWts{1}, 0);     
gH_Params.PWts = PWts;

% % create Vandermonde matrix
% nxi = intParams.nrPoints ;
% nn = intParams.d ;
% V = [ones(nxi,1), repmat(pts,1,nn).^repmat(1:nn,nxi,1)] ;
% [rV,cV] = size(V) ;

application_matrix = vector_partial_application(1, p, intParams.mon_basis);
probData.A = sparse([-eye(U), intParams.mon_to_P0 * application_matrix*intParams.P0_to_mon]);
eval_rhs = msspoly(c);
eval_rhs_vec = msspoly_to_vector(eval_rhs, intParams.mon_basis);
probData.c = [zeros(U, 1); -intParams.mon_to_P0 * eval_rhs_vec]; 
probData.b = -intParams.w;

% % create A,b,c matrices to define conic problem
% % probData.A = sparse([eye(U),[V,zeros(nxi,U-cV)]']) ;
% probData.A = sparse([eye(U),ones(U,1),zeros(U,U-1)]) ;
% probData.b = intParams.w ;
% probData.c = -[zeros(U,1); c; -10000*ones(U-1,1)] ;

% make initial dual iterate
y0 = ones(size(probData.c)) ;
[~, g0, ~, ~] = alfonso_grad_and_hess(y0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*y0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
y0 = repmat(sqrt(rP*rD),2*U,1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, y0, @alfonso_grad_and_hess, gH_Params, opts);

%% PLOTTING RESULTS
figure(1) ; cla ; hold on ;
xlim([-1, 1])
ylim([-1, 2])

tvec = linspace(-1,1,500) ;

% spotless output
fvals = msubs(fmosek,t,tvec) ;
plot(tvec,fvals,'LineWidth',1.5)

results.y
% alfonso output:
%yvals = -results.y ; % these are the values of the polynomial f at the points "pts"
mon = (intParams.P0_to_mon * results.y)' * intParams.mon_basis.monomials
monvals = msubs(mon,intParams.mon_basis.variables,tvec) ;
plot(tvec, monvals,'--','LineWidth',2.5)

% %% ALFONSO GRADIENT FUNCTIONS
% function [in, g, H, L] = gH_SOS_at_a_point(x, params)
%     in = 1;
%     U = params.U;
%     P = params.P;
%     PWts = params.PWts;
%     p = params.p;
%     c = params.c;
%     
% % first constraint: f >= 0
%     [inPoly,gPoly,HPoly] = gH_SOSWt(x(1:U),P) ;
% 
%     if inPoly == 1
%         % for the weight 1-t_j^2
%         [inPolyWt, gPolyWt, HPolyWt] = gH_SOSWt(x(1:U),PWts);
%         inPoly  = inPoly & inPolyWt;
%         if inPoly == 1
%             gPoly   = gPoly+gPolyWt;
%             HPoly   = HPoly+HPolyWt;
%         else
%             gPoly   = NaN;
%             HPoly   = NaN;
%         end
%     end
%     
% % second constraint: f(p) >= 0
%     [~, gp, Hp] = gH_SOSWt(x(end),c) ;
%     g = [gPoly;gp] ;
%     H = blkdiag(HPoly,Hp) ;
%     
%     % check that H is PSD
%     [L,err] = chol(H,'lower') ;
%     if err ~= 0
%         g = NaN ;
%         H = NaN ;
%         L = NaN ;
%     end
% end
% 
% 
% function [in, g, H] = gH_SOSWt(x, P)
% % This method computes the gradient and Hessian of a barrier function term
% % corresponding to a single weight and single approximated polynomial for
% % the problem of polynomial envelopes.
% % --------------------------------------------------------------------------
% % USAGE of "gH_SOSWt"
% % [in, g, H] = gH_SOSWt(x, P)
% % --------------------------------------------------------------------------
% % INPUT
% % x:    subvector of the primal iterate corresponding to a single
% %       approximated polynomial
% % P:    evaluations of "weighted" basis polynomials at the
% %       interpolation points
% %
% % OUTPUT
% % in:   0 if P'*diag(x)*P is not positive definite. 1 if P'*diag(x)*P is
% %       positive definite.
% % g:    gradient of the barrier function term corresponding to a single
% %       weight and single approximated polynomial at x
% % H:    Hessian of the barrier function term corresponding to a single
% %       weight and single approximated polynomial at x
% % --------------------------------------------------------------------------
% % EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% % None.
% % -------------------------------------------------------------------------
% 
%     Y = P'*diag(x)*P;
% 
%     if ~issymmetric(Y)
%         Y = (Y+Y')/2;
%     end
%     
%     [L, err] = chol(Y, 'lower');
%     if err > 0
%         in = 0;
%         g = NaN;
%         H = NaN;
%     else
%         in = 1;
%         V = L\P';
%         VtV = V'*V;
% 
%         g = -diag(VtV);
%         H = VtV.^2;
%     end
% end