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
function example_program()

variables = msspoly('x');
fprintf("
f=variables(1)^3+3;
degree = 10;
S_bounds = [-1, 2];


hS = -(variables-S_bounds(:, 1)).*(variables-S_bounds(:, 2));


dl=boxMoments(variables, S_bounds(:, 1), S_bounds(:, 2));

prog = AlfonsoSOSProgFekete;
prog.with_indeterminate(variables);

[w, wcoeff, w_monomial] = prog.new_free_poly(variables, degree) ;
bad_operator = diff(w, variables) * f;

prog.sos_on_K(bad_operator-1, variables, S_bounds, degree+ceil(msspoly_degree(f)/2) * 2);
prog.sos_on_K(w, variables, S_bounds, degree);
obj = dl(wmonom)' * wcoeff ;

sol_alfonsless = prog.minimize(obj);

l2_dist = l2_dist_on_box(fspotless, sol_alfonsless.polys(1), S_bounds(:, 1), S_bounds(:, 2), variables)



% close all
% figure('Name', 'dw/dt * f') ; cla ; hold on ;
% 
% tvec = linspace(-1,2,500) ;
% 
% const2_vals = dmsubs(constraint2, variables, tvec);
% const1_vals = dmsubs(constraint1, variables, tvec);
% 
% 
% plot(tvec, const2_vals, '-k')
% plot(tvec, const1_vals, '--r')
% 
% 
% figure('Name', 'w') ; cla ; hold on ;
% %xlim([-1, 2])
% %ylim([-1, 2])
% 
% tvec = linspace(-1,2,500) ;
% 
% % spotless output
% fvals = dmsubs(fspotless,variables,tvec) ;
% plot(tvec,fvals,'--', 'LineWidth',1.5)
% 
% % alfonso output:
% %yvals = -results.y ; % these are the values of the polynomial f at the points "pts"
% 
% monvals = dmsubs(falfonso,variables,tvec);
% plot(tvec,monvals,'LineWidth',1);
% 
% end
