 %% problem description
% This does not work correctly. Something about how the A matrix is defined
% is broken.
%
%   min   int f(x,y) dxdy
%   f,g   X*Y
%
%   s.t.    f     >=  g     on X*Y
%           g     >=  0     on X*Y
%         g(0,.)  >=  1     on Y
%         f(.,0)  >=  1     on X
%
%% user parameters
degree = 12 ;

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
[prog,g,gcoeff] = prog.newFreePoly(mon) ;

% create SOS constraints
gcon = msubs(g,x,0) ;
fcon = msubs(f,y,0) ;
prog = sosOnK(prog, f - g, [x;y], [hX;hY], degree) ;
prog = sosOnK(prog, g, [x;y], [hX;hY], degree) ;
prog = sosOnK(prog, gcon - 1 , y, hY, degree) ;
prog = sosOnK(prog, fcon - 1 , x, hX, degree) ;

% create cost function
int_XY = boxMoments([x;y],[-1;-1],[1;1]) ;
obj = int_XY(mon)'*fcoeff ;

% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;
sol = prog.minimize(obj, @spot_mosek, options);

% recover f and g
fspotless = sol.eval(fcoeff)'*mon ;
gspotless = sol.eval(gcoeff)'*mon ;

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
numPolys = 4 ;  % number of constraints
gH_Params.n = n;
gH_Params.d = d;
gH_Params.U = U;
gH_Params.numPolys = numPolys;

gH_Params.L     = L;
gH_Params.LWts  = LWts;
nu              = numPolys*(L+sum(LWts));
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

% find the interpolation points where x = 0; we fill in the cone problem's
% "c" matrix with the column vectors produced here
plog = abs(pts(:,1) - 0) < 1e-6 ;
xconval = -10*ones(U,1) ;
xconval(plog) = 1 ;

% find the interpolation points where y = 0
plog = abs(pts(:,2) - 0) < 1e-6 ;
yconval = -10*ones(U,1) ;
yconval(plog) = 1 ;

% create A and c:
% constraint 1: f - g >= 0
A1 = sparse([eye(U), -eye(U)]) ;
c1 = zeros(U,1) ;

% constraint 2: g >= 0
A2 = sparse([zeros(U), eye(U)]) ;
c2 = zeros(U,1) ;

% constraint 3: g(0,.) >= 1
A3 = sparse([zeros(U), eye(U)]) ;
c3 = yconval ;

% constraint 4: f(.,0) >= 1
A4 = sparse([eye(U), zeros(U)]) ;
c4 = xconval ;

% create A,b,c matrices to define conic problem
probData.A = [A1 ; A2 ; A3 ; A4]' ; % <-- this does not work!
probData.b = repmat(intParams.w,2,1) ;
probData.c = [c1;c2;c3;c4] ;

% make initial primal iterate
x0 = ones(size(probData.c)) ;
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),4*U,1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

%% RECOVER ALFONSO SOLUTION AS MSSPOLY
% recover solution
zvals = -results.y(1:end/2) ; % the 1st half is for f, the 2nd for g
data = msubs(mon,[x;y],pts') ;
coeff = data'\zvals ;
falfonso = coeff'*mon ;

%% PLOTTING POLYNOMIAL (RECOVERED) OUTPUT
figure(1) ; cla ;

% spotless setup
[X,Y] = makeContourAxes() ;
XY = makeContourAxes() ;

% spotless output
F = reshape(full(msubs(fspotless,[x;y],XY)),100,100) ;
subplot(1,2,1)
surf(X,Y,F)

% spotless output
F = reshape(full(msubs(falfonso,[x;y],XY)),100,100) ;
subplot(1,2,2) ;
surf(X,Y,F)

title('The left and right subplots should be the same :(')