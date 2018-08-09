%% problem description
%
%   min     int  f(x,y) dxdy
%    f      X*Y
%
%   s.t.    f(x,y)  >= -1   for all (x,y) \in X*Y
%           f(c,.)  >= +1   for all y \in Y
%
%% user parameters 
degree = 12 ;
c = 0 ;

%% SPOTLESS PROBLEM
% variables
x = msspoly('x',1) ;
y = msspoly('y',1) ;

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(x) ;
prog = prog.withIndeterminate(y) ;

% bounds on variables
hX = (x + 1).*(1 - x) ;
hY = (y + 1).*(1 - y) ;

% create decision variable
mon = monomials([x;y],0:degree) ;
[prog,f,fcoeff] = prog.newFreePoly(mon) ;

% create SOS constraints
xcon = msubs(f,x,c) ;
prog = sosOnK(prog,f+1,[x;y],[hX;hY],degree) ;
prog = sosOnK(prog,xcon-1,y,hY,degree) ;

% create cost function
int_XY = boxMoments([x;y],[-1;-1],[1;1]) ;
obj = int_XY(mon)'*fcoeff ;
int_XY
assert 0;

% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;
sol = prog.minimize(obj, @spot_mosek, options);

% recover f
fspotless = sol.eval(fcoeff)'*mon ;

%% ALFONSO PROBLEM
% create interpolation points
intParams = FeketeCube(2,degree/2) ;
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
PWts = cell(n,1);
for j = 1:n
    PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
    [PWts{j}, ~]    = qr(PWts{j}, 0);
    % associated positive semidefinite cone constraints: 
    % PWts{j}'*diag(x_1)*PWts{j} >= 0,
    % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
end
gH_Params.PWts = PWts;

% find the points where x = c
plog = abs(pts(:,1) - c) < 1e-6 ;
cpval = -10*ones(U,1) ;
cpval(plog) = 1 ;

% create A,b,c matrices to define conic problem
probData.A = sparse(repmat(eye(U),1,2)) ;
probData.b = intParams.w ;
probData.c = -[-ones(U,1) ; cpval] ;

% make initial primal iterate
x0 = ones(2*U,1) ;
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),2*U,1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

%% RECOVER ALFONSO SOLUTION AS MSSPOLY
% recover solution
zvals = -results.y ;
data = msubs(mon,[x;y],pts') ;
coeff = data'\zvals ;
falfonso = coeff'*mon ;

%% COMPARE SOLUTIONS
fvals_spotless = msubs(fspotless,[x;y],pts') ;

% fvals_alfonso  = -results.y ; % the direct output from alfonso, should
%                               % have nearly identical values to the line
%                               % below

fvals_alfonso =  msubs(falfonso,[x;y],pts') ;

max(abs(fvals_spotless(:) - fvals_alfonso(:)))

%% PLOTTING POLYNOMIAL (RECOVERED) OUTPUT
% spotless setup
[X,Y] = makeContourAxes() ;
XY = makeContourAxes() ;

% alfonso setup
xvals = pts(:,1) ;
yvals = pts(:,2) ;
tri = delaunay(xvals,yvals);

figure(1) ; cla ; hold on ;

% spotless output
F1 = reshape(full(msubs(fspotless,[x;y],XY)),100,100) ;
% surf(X,Y,F)

% alfonso output
F2 = reshape(full(msubs(falfonso,[x;y],XY)),100,100) ;

% plot error between spotless and alfonso on X*Y
surf(X,Y,F1-F2)

%% PLOTTING COARSE OUTPUT
figure(2) ; cla ; hold on

% spotless output
zvals = full(fvals_spotless) ;
trisurf(tri, xvals,yvals,zvals);

% alfonso output
zvals = -results.y ;
trisurf(tri, xvals,yvals,zvals);
