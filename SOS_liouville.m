%% problem description
%
%   min     int  f(x,y) dxdy
%    f      X*Y
%
%   s.t.    f(x,y)  >= -1   for all (x,y) \in X*Y
%           f(c,.)  >= +1   for all y \in Y
%
%% user parameters 
degree = 8 ;
c = 0 ;
t = msspoly('t');
x = msspoly('x');
variables = [t;x];
f = variables(2)^2;
T = 1;
intParams_w = FeketeCube(size(x, 1), degree/2,x) ;

intParams_v_f = FeketeCube(size(variables, 1), degree/2 + ceil(msspoly_degree(f)/2), variables); 
intParams_arr = [intParams_v_f, intParams_w, intParams_w, intParams_w];

intParams_v = FeketeCube(size(variables, 1), degree/2, variables);
% %% SPOTLESS PROBLEM
% % variables
% x = intParams.mon_basis.variables(1) ;
% y = intParams.mon_basis.variables(2) ;

% % create program
% prog = spotsosprog ;
% prog = prog.withIndeterminate(x) ;
% prog = prog.withIndeterminate(y) ;

% % bounds on variables
% hX = (x + 1).*(1 - x) ;
% hY = (y + 1).*(1 - y) ;

% % create decision variable
% mon = monomials([x;y],0:degree) ;
% [prog,f,fcoeff] = prog.newFreePoly(mon) ;

% % create SOS constraints
% xcon = msubs(f,x,c) ;
% prog = sosOnK(prog,f+1,[x;y],[hX;hY],degree) ;
% prog = sosOnK(prog,xcon-1,y,hY,degree) ;
% prog = sosOnK(prog, f, [x;y], [hX; hY], degree)

% % create cost function
% int_XY = boxMoments([x;y],[-1;-1],[1;1]) ;
% obj = int_XY(mon)'*fcoeff ;

% % run solver
% options = spot_sdp_default_options();
% options.verbose = 1 ;
% sol = prog.minimize(obj, @spot_mosek, options);

% % recover f
% fspotless = sol.eval(fcoeff)'*mon ;

%% ALFONSO PROBLEM
% create interpolation points

numPolys = 4;
gH_Params.numPolys = numPolys;
gH_Params.n_arr = zeros(numPolys, 1);
gH_Params.d_arr = zeros(numPolys, 1);
gH_Params.U_arr = zeros(numPolys, 1);
gH_Params.L_arr = zeros(numPolys, 1);
gH_Params.LWts_cell = cell(numPolys, 1);
gH_Params.P_cell = cell(numPolys, 1);
gH_Params.PWts_cell = cell(numPolys, 1);
gH_Params.bnu = 0;

for i=1:numPolys
    intParams_i = intParams_arr(i);
    n   = intParams_i.n;
    d   = intParams_i.d;
    U   = intParams_i.U;
    L   = intParams_i.L;
    P   = intParams_i.P0;
    pts = intParams_i.pts;

    % dimension of the weight polynomial space (should be dimension of d)
    LWts = repmat(nchoosek(n+d-1,n),n,1);

    % parameter object for cone gradient/hessian function
    gH_Params.n_arr(i)      = n;
    gH_Params.d_arr(i)      = d;
    gH_Params.U_arr(i)      = U;
    gH_Params.L_arr(i)      = L;
    gH_Params.LWts_cell{i}  = LWts;
    nu                      = L+sum(LWts) ;
    gH_Params.bnu           = gH_Params.bnu + nu;
    gH_Params.P_cell{i}     = P;


    % create polynomial hT (g in the alfonso paper) to define space T = [-1,1]^2
    wtVals = 1-pts.^2;
    PWts = cell(n,1);
    for j = 1:n
        PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
        [PWts{j}, ~]    = qr(PWts{j}, 0);
        % associated positive semidefinite cone constraints: 
        % PWts{j}'*diag(x_1)*PWts{j} >= 0,
        % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
    end
    gH_Params.PWts_cell{i} = PWts;
end
gH_Params.bnu = gH_Params.bnu + 1;

v_to_v_f = monomial_to_monomial(intParams_v.mon_basis, intParams_v_f.mon_basis);
time_der_v = vector_derivative(1, intParams_v.mon_basis);
mult_mat = vector_poly_multiply(f, intParams_v.mon_basis, intParams_v_f.mon_basis);
space_der_v = vector_derivative(2, intParams_v.mon_basis);

A1 = intParams_v_f.mon_to_P0 * (v_to_v_f * time_der_v + mult_mat * space_der_v) * intParams_v.P0_to_mon;
A2 = zeros(intParams_v_f.U, intParams_w.U);

v_to_w = monomial_to_monomial(intParams_v.mon_basis, intParams_w.mon_basis);
apply_final_time = vector_partial_application(1, T, intParams_v.mon_basis);

B1 = intParams_w.mon_to_P0 * v_to_w * apply_final_time * intParams_v.P0_to_mon;
B2 = zeros(intParams_w.U);

apply_init_time = vector_partial_application(1, -1, intParams_v.mon_basis);

C1 = intParams_w.mon_to_P0 * v_to_w * apply_init_time * intParams_v.P0_to_mon;
C2 = eye(intParams_w.U);

D1 = zeros(intParams_w.U, intParams_v.U);
D2 = eye(intParams_w.U);

const1 = zeros(intParams_v_f.U, 1);
const2 = zeros(intParams_w.U, 1);
const3 = ones(intParams_w.U, 1);
const4 = zeros(intParams_w.U, 1);

int1 = zeros(intParams_v.U, 1);
int2 = intParams_w.w;

A = -[A1, A2;
     B1, B2;
     C1, C2;
     D1, D2]';
c = -[const1;
      const2;
      const3;
      const4];
b = -[int1;
      int2];

% find the points where x = c
% plog = abs(pts(:,1) - c) < 1e-6 ;
% cpval = -10*ones(U,1) ;
% cpval(plog) = 1 ;

% create A,b,c matrices to define conic problem
% probData.A = sparse(repmat(eye(U),1,2)) ;
% probData.b = intParams.w ;
% probData.c = -[-ones(U,1) ; cpval] ;
% make initial primal iterate
x0 = ones(gH_Params.numPolys*U, 1);
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),gH_Params.numPolys*U,1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

%% RECOVER ALFONSO SOLUTION AS MSSPOLY
% recover solution
%zvals = -results.y ;
%data = msubs(mon,[x;y],pts') ;
%coeff = data'\zvals ;
%falfonso = coeff'*mon ;

falfonso = (intParams.P0_to_mon * results.x(1:U, :))' * intParams.mon_basis.monomials;

%% COMPARE SOLUTIONS
fvals_spotless = msubs(fspotless,[x;y],pts') ;

% fvals_alfonso  = -results.y ; % the direct output from alfonso, should
%                               % have nearly identical values to the line
%                               % below

fvals_alfonso =  msubs(falfonso,[x;y],pts') ;

%max(abs(fvals_spotless(:) - fvals_alfonso(:)))

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
