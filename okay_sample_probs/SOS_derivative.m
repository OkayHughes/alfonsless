%% problem description
%
%   min     int  f(x,y) dxdy
%    f      X*Y
%
%   s.t.    f(x,y)  >= -1   for all (x,y) \in X*Y
%           f(c,.)  >= +1   for all y \in Y
%
%% user parameters 


%% ALFONSO PROBLEM
% create interpolation points
function SOS_derivative()

variables = msspoly('x');
f=variables(1)^3+3;
degree = 10;
S_bounds = [-1, 2];
intParams_w = FeketeCube(size(variables, 1), degree/2, variables);
intParams_w = scale_fekete_cube(intParams_w, S_bounds);


intParams_w_f = FeketeCube(size(variables, 1), degree/2 + ceil(msspoly_degree(f)/2), variables); 
intParams_w_f = scale_fekete_cube(intParams_w_f, S_bounds);

intParams_arr = [intParams_w, intParams_w_f];

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
    gH_Params.bnu           = gH_Params.bnu + nu;
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



% create A,b,c matrices to define conic problem
% probData.A = sparse(repmat(eye(U),1,2)) ;
% probData.b = intParams.w ;
% probData.c = -[-ones(U,1) ; cpval] ;

A1 = eye(intParams_w.U);

mult_mat = vector_poly_multiply(f, intParams_w.mon_basis, intParams_w_f.mon_basis);
time_der = vector_derivative(1, intParams_w.mon_basis);

B1 = intParams_w_f.mon_to_P0 * mult_mat * time_der * intParams_w.P0_to_mon;

const1 = zeros(intParams_w.U, 1);
const2 = ones(intParams_w_f.U, 1);

int1 = intParams_w.w;

probData.A = -[A1;
              B1]';

probData.c = -[const1;  
              const2];
probData.b = -[int1];
% make initial primal iterate
x0 = ones(sum(gH_Params.U_arr), 1);
[~, g0, ~, ~] = alfonso_grad_and_hess(x0, gH_Params);
rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
rD = max((1+abs(g0))./(1+abs(probData.c)));
x0 = repmat(sqrt(rP*rD),sum(gH_Params.U_arr),1);  

% run alfonso
opts.optimTol = 1e-6 ;
results = alfonso(probData, x0, @alfonso_grad_and_hess, gH_Params, opts);

falfonso = (intParams_w.P0_to_mon * results.y)' * intParams_w.mon_basis.monomials;

%run spotless problem
hS = -(variables-S_bounds(:, 1)).*(variables-S_bounds(:, 2));
dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));
prog = spotsosprog;
prog = prog.withIndeterminate(variables);
wmonom = monomials(variables, 0:degree) ;
[prog, w, wcoeff] = prog.newFreePoly(wmonom) ;
bad_operator = diff(w, variables) * f;
prog = sosOnK(prog, bad_operator-1, variables, hS, degree) ;
prog = sosOnK(prog, w, variables, hS, degree);
obj = dl(wmonom)' * wcoeff ;

options = spot_sdp_default_options();
options.verbose = 1;
sol = prog.minimize(obj, @spot_mosek, options);

fspotless = sol.eval(w) ;





figure(1) ; cla ; hold on ;
%xlim([-1, 2])
%ylim([-1, 2])

tvec = linspace(-1,2,500) ;

% spotless output
fvals = dmsubs(fspotless,variables,tvec) ;
plot(tvec,fvals,'--', 'LineWidth',1.5)

% alfonso output:
%yvals = -results.y ; % these are the values of the polynomial f at the points "pts"

monvals = dmsubs(falfonso,variables,tvec);
plot(tvec,monvals,'LineWidth',1);

end
