%similar to SOS_on_a_subset_okay, except that it uses the new formulation 
%of the hessian, and s_2 (the second slack variable) is an element of the
% univariate WSOS cone

%% problem description
%
%   min     int  f(x,y) dxdy
%    f      X*Y
%
%   s.t.    f(x,y)  >= -1   for all (x,y) \in X*Y
%           f(c,.)  >= +1   for all y \in A \subset Y
%
%% user parameters 
function SOS_on_an_affine_2d_okay()

degree = 12;
c = 0 ;

sp_bounds = [-1, 1; -1, 1];
A_bounds = [-0.2, 0.2;
            -0.2, 0.2];

intParams = FeketeCube(2,degree/2) ;
intParams = scale_fekete_cube(intParams, sp_bounds);
sintParams = FeketeCube(2, degree/2, intParams.mon_basis.variables);
sintParams = scale_fekete_cube(sintParams, A_bounds);

intParams_arr = [intParams, sintParams];
% %% SPOTLESS PROBLEM
% variables
x = intParams.mon_basis.variables(1) ;
y = intParams.mon_basis.variables(2) ;

% create program
prog = spotsosprog ;
prog = prog.withIndeterminate(x) ;
prog = prog.withIndeterminate(y) ;

% bounds on variables
vars = [x;y]
hs = -(vars-sp_bounds(:, 1)).*(vars-sp_bounds(:, 2));
hX = hs(1);
hY = hs(2);
ha = -(vars-A_bounds(:, 1)).*(vars-A_bounds(:, 2));


% create decision variable
mon = monomials([x;y],0:degree) ;
[prog,f,fcoeff] = prog.newFreePoly(mon) ;

% create SOS constraints
%prog = sosOnK(prog,f+1,[x;y],[hX;hY],degree) ;
prog = sosOnK(prog,f-1,vars,ha,degree) ;
prog = sosOnK(prog, f+1, [x;y], [hX; hY], degree)

% create cost function
int_XY = boxMoments([x;y],sp_bounds(:, 1),sp_bounds(:, 2)) ;
obj = int_XY(mon)'*fcoeff ;

% run solver
options = spot_sdp_default_options();
options.verbose = 1 ;
sol = prog.minimize(obj, @spot_mosek, options);

% recover f
fspotless = sol.eval(fcoeff)'*mon ;

%% ALFONSO PROBLEM
% create interpolation points
numPolys = 2;
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
    bounds = intParams_i.bounds;
    lb = bounds(:, 1);
    ub = bounds(:, 2);

    % dimension of the weight polynomial space (should be dimension of d)
    LWts = repmat(nchoosek(n+d-1,n),n,1);

    % parameter object for cone gradient/hessian function
    gH_Params.n_arr(i)      = n;
    gH_Params.d_arr(i)      = d;
    gH_Params.U_arr(i)      = U;
    gH_Params.L_arr(i)      = L;
    gH_Params.LWts_cell{i}  = LWts;
    nu                      = L+sum(LWts) ;
    gH_Params.bnu	        = gH_Params.bnu + nu;
    gH_Params.P_cell{i}     = P;


    % create polynomial hT (g in the alfonso paper) to define space T = [-1,1]^2
    wtVals  = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);
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

% find the points where x = c
% plog = abs(pts(:,1) - c) < 1e-6 ;
% cpval = -10*ones(U,1) ;
% cpval(plog) = 1 ;

% create A,b,c matrices to define conic problem
% probData.A = sparse(repmat(eye(U),1,2)) ;
% probData.b = intParams.w ;
% probData.c = -[-ones(U,1) ; cpval] ;
application_matrix = vector_partial_application(1, c, intParams.mon_basis);
l2s = monomial_to_monomial(intParams.mon_basis, sintParams.mon_basis);

probData.A = sparse([-eye(gH_Params.U_arr(1)); -sintParams.mon_to_P0 * l2s *  intParams.P0_to_mon]');
mon_con1 = msspoly(-1);
mon_con1_vec = intParams.mon_to_P0 * msspoly_to_vector(mon_con1, intParams.mon_basis);
mon_con2 = msspoly(1);
mon_con2_vec = sintParams.mon_to_P0 * msspoly_to_vector(mon_con2, sintParams.mon_basis);
probData.c = [-mon_con1_vec; -mon_con2_vec];
probData.b = -[intParams.w];
% make initial primal iterate
x0 = ones(sum(gH_Params.U_arr), 1);
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD), sum(gH_Params.U_arr),1);  

% run alfonso
opts.optimTol = 1e-6 ;
%opts.verbose = 0;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

spot_poly_vec = msspoly_to_vector(fspotless, intParams.mon_basis);
alfonso_poly_vec = intParams.P0_to_mon * results.y;
falfonso = alfonso_poly_vec' * intParams.mon_basis.monomials;

[grid_x, grid_y] = meshgrid(linspace(-1, 1, 50), linspace(-1, 1, 50));
flat = [grid_x(:), grid_y(:)];
alf_vals = reshape(dmsubs(falfonso, intParams.mon_basis.variables, flat')', 50, 50);
spot_vals = reshape(dmsubs(fspotless, intParams.mon_basis.variables, flat')', 50, 50);


dist_func = l2_dist_on_box(falfonso, fspotless, sp_bounds(:, 1), sp_bounds(:, 2), vars)

COA(:,:,1) = zeros(50); % red
COA(:,:,2) = zeros(50); % green
COA(:,:,3) = zeros(50); % blue

COS(:,:,1) = ones(50); % red
COS(:,:,2) = ones(50); % green
COS(:,:,3) = ones(50); % blue




figure ; cla ; hold on;

surf(grid_x, grid_y, alf_vals, COA);
surf(grid_x, grid_y, spot_vals, COS);

end

%% RECOVER ALFONSO SOLUTION AS MSSPOLY
% % recover solution
% %zvals = -results.y ;
% %data = msubs(mon,[x;y],pts') ;
% %coeff = data'\zvals ;
% %falfonso = coeff'*mon ;

% falfonso = (intParams.P0_to_mon * results.x(1:U, :))' * intParams.mon_basis.monomials;

% %% COMPARE SOLUTIONS
% fvals_spotless = msubs(fspotless,[x;y],pts') ;

% % fvals_alfonso  = -results.y ; % the direct output from alfonso, should
% %                               % have nearly identical values to the line
% %                               % below

% fvals_alfonso =  msubs(falfonso,[x;y],pts') ;


%max(abs(fvals_spotless(:) - fvals_alfonso(:)))

% %% PLOTTING POLYNOMIAL (RECOVERED) OUTPUT
% % spotless setup
% [X,Y] = makeContourAxes() ;
% XY = makeContourAxes() ;

% % alfonso setup
% xvals = pts(:,1) ;
% yvals = pts(:,2) ;
% tri = delaunay(xvals,yvals);

% figure(1) ; cla ; hold on ;

% % spotless output
% F1 = reshape(full(msubs(fspotless,[x;y],XY)),100,100) ;
% % surf(X,Y,F)

% % alfonso output
% F2 = reshape(full(msubs(falfonso,[x;y],XY)),100,100) ;

% % plot error between spotless and alfonso on X*Y
% surf(X,Y,F1-F2)

% %% PLOTTING COARSE OUTPUT
% figure(2) ; cla ; hold on

% % spotless output
% zvals = full(fvals_spotless) ;
% trisurf(tri, xvals,yvals,zvals);

% % alfonso output
% zvals = -results.y ;
% trisurf(tri, xvals,yvals,zvals);
