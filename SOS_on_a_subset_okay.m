%% problem description
%
%   min     int  f(x,y) dxdy
%    f      X*Y
%
%   s.t.    f(x,y)  >= -1   for all (x,y) \in X*Y
%           f(c,1)  >= +1   for all y \in Y
%
%% user parameters 
degree = 12 ;
c = 0 ;
intParams1 = FeketeCube(2,degree);
intParams2 = FeketeCube(1, degree, intParams1.mon_basis.variables(2));
intParams_arr = [intParams1, intParams2];

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
numPolys = 2;
gH_Params.numPolys = numPolys;
L = intParams_arr(1).L;
gH_Params.L = intParams_arr(1).L;
gH_Params.n_arr = zeros(numPolys, 1);
gH_Params.d_arr = zeros(numPolys, 1);
gH_Params.U_arr = zeros(numPolys, 1);
gH_Params.LWts_arr = cell(numPolys, 1);
gH_Params.P_cell = cell(numPolys, 1);
gH_Params.Pwts_cell = cell(numPolys, 1);

for i=1:numPolys
    intParams = intParams_arr(i);
    n   = intParams.n;
    d   = intParams.d;
    U   = intParams.U;
    P   = intParams.P0;
    pts = intParams.pts;

    % parameter object for cone gradient/hessian function
    gH_Params.n_arr(i) = n;
    gH_Params.d_arr(i) = d;
    gH_Params.U_arr(i) = U;


    gH_Params.P_cell{i} = P;
    % dimension of the weight polynomial space (should be dimension of d)
    LWts_i = nchoosek(n+d-1,n) * ones(n, 1);
    gH_Params.LWts_arr{i} = LWts_i; %todo change

end

for i=1:numPolys
    intParams = intParams_arr(i);
    n   = intParams.n;
    d   = intParams.d;
    U   = intParams.U;
    P   = intParams.P0;
    pts = intParams.pts;
    LWts = gH_Params.LWts_arr{j}
    i = i
    % create polynomial hT (g in the alfonso paper) to define space T = [-1,1]^2
    wtVals = 1-pts.^2;
    PWts = cell(n,1);
    for j = 1:n
        j = j
        size(P)
        gH_Params.LWts_arr(j)
        PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
        [PWts{j}, ~]    = qr(PWts{j}, 0);
        % associated positive semidefinite cone constraints: 
        % PWts{j}'*diag(x_1)*PWts{j} >= 0,
        % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
    end
    gH_Params.PWts_cell{i} = PWts;
end

gH_Params.bnu   = (sum(sum(cell2mat(gH_Params.LWts_arr))) + gH_Params.numPolys * gH_Params.L) + 1;
% find the points where x = c
% plog = abs(pts(:,1) - c) < 1e-6 ;
% cpval = -10*ones(U,1) ;
% cpval(plog) = 1 ;

% create A,b,c matrices to define conic problem
% probData.A = sparse(repmat(eye(U),1,2)) ;
% probData.b = intParams.w ;
% probData.c = -[-ones(U,1) ; cpval] ;
application_matrix = vector_partial_application(1, c, intParams_arr(1).mon_basis);
mon_change_basis = monomial_to_monomial(intParams_arr(1).mon_basis, intParams_arr(2).mon_basis);
size(-eye(intParams_arr(1).U))
size(intParams_arr(2).mon_to_P0 * mon_change_basis * application_matrix * intParams_arr(1).P0_to_mon)
probData.A = sparse([-eye(intParams_arr(1).U); ...
                     (-intParams_arr(2).mon_to_P0 * mon_change_basis * application_matrix * intParams_arr(1).P0_to_mon)])';
mon_con1 = msspoly(-1);
mon_con1_vec = intParams_arr(1).mon_to_P0 * msspoly_to_vector(mon_con1, intParams_arr(1).mon_basis);
mon_con2 = msspoly(1);
mon_con2_vec = intParams_arr(2).mon_to_P0 * msspoly_to_vector(mon_con2, intParams_arr(2).mon_basis);
probData.c = [-mon_con1_vec; -mon_con2_vec];
probData.b = -[intParams_arr(1).w];
% make initial primal iterate
x0 = ones(sum(gH_Params.U_arr), 1);
size(probData.A)
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
size(probData.A*x0)
size(probData.b)
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),gH_Params.numPolys*U,1);  

% run alfonso
opts.optimTol = 1e-6;
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
