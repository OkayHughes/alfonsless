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

intParams_w = AffineFeketeCube(variables, degree, S_bounds);


intParams_w_f = AffineFeketeCube(variables, degree + ceil(msspoly_degree(f)/2)*2, S_bounds);

intParams_arr = [intParams_w, intParams_w_f];

gH_Params = gen_grad_params(intParams_arr);



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


constraint1 = (mult_mat * time_der * intParams_w.P0_to_mon * results.y)' * intParams_w_f.mon_basis.monomials;
constraint2 = diff(falfonso, variables) * f;

close all
figure('Name', 'dw/dt * f') ; cla ; hold on ;

tvec = linspace(-1,2,500) ;

const2_vals = dmsubs(constraint2, variables, tvec);
const1_vals = dmsubs(constraint1, variables, tvec);


plot(tvec, const2_vals, '-k')
plot(tvec, const1_vals, '--r')


figure('Name', 'w') ; cla ; hold on ;
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
