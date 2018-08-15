clear
close all

% options
degree = 8 ;
solver = 'mosek' ;
sostype = 'sos' ;
visualize = 1 ;
plotting_color = [1 0 0] ;

%% set up advection vars
t = msspoly('t', 1);
x = msspoly('x', 1);
k = msspoly('k', 1) ;

%% set up spaces
T = 1 ;

X_range = repmat([-1, 1], size(x)) ;

X0_range = repmat([-0.2, 0.2], size(x)) ;

K_range = repmat([-1,1], size(k)) ;

%% dynamics
% system dynamics
f = 0.2*k ;
g = 0.5*k^2 ;
%g = 0*k;
% g = 4.4*(t-0.5)^2  ;
% g = -4*(t-0.5)^2 + 1.0  ;

%% create domain
hX = (x - X_range(:,1)).*(X_range(:,2) - x) ;

hX0 = (x - X0_range(:,1)).*(X0_range(:,2) - x) ;
% hX0 = 1 - (x(1)/X0_range(1,1))^2 - (x(2)/X0_range(1,2))^2 ;

hK = (k - K_range(:,1)).*(K_range(:,2) - k) ;

% pin x0 to 0
%hDelta = [1 0] ;

%% create cost function
int_XK = boxMoments([x;k], [X_range(:,1); K_range(:,1)], [X_range(:,2); K_range(:,2)]) ;

%% setup the problem structure
prob.t = t ;
prob.z = x ;
prob.k = k ;
prob.f = f ;
prob.g = g ;
prob.hZ = hX ;
prob.hZ0 = hX0 ;
prob.hK = hK ;
%prob.hDelta = hDelta ;
prob.int_ZK = int_XK ;
prob.degree = degree ;
prob.domain_size = 1 ;
prob.solver = solver ;
prob.run_solver = 1 ;

%% solve the problem
out = findFRS(prob) ;



prob_a.degree = degree;
prob_a.t = t;
prob_a.x = x;
prob_a.k = k;
prob_a.f = f;
prob_a.g = g;
prob_a.X_bounds = X_range;
prob_a.X0_bounds = X0_range;
prob_a.K_bounds = K_range;
prob_a.T = T;

out_a = find_FRS_alfonso(prob_a);

Y_bounds = [X_range; K_range];
l2_dist = l2_dist_on_box(out.w, out_a.w, Y_bounds(:, 1), Y_bounds(:, 2), [x;k])
            

%% set up variables for visualization
N = 100 ;

xvec = linspace(X_range(:,1),X_range(:,2),N) ;
kvec = linspace(K_range(:,1),K_range(:,2),N) ;

[X,K] = meshgrid(xvec, kvec) ;
XK = [X(:) K(:)]' ;


%% output visualization
if visualize
%     v = out.v ;
% 
%     N = 100 ;
% 
%     close all
%     figure(1)
% 
%     % subplot(2,1,1)
%     hold on
%     tic
%     for tidx = 0:0.25:1
%         v0 = msubs(v,t,tidx) ;
% 
%         F = msubs(v0, [x;k], XK) ;
%         F = full(F) ;
%         F = reshape(F,N,N) ;
% 
%         contour(xvec,kvec,F,[0 0],'Color',[tidx, 0.5, 1-tidx])
%     end
%     axis equal
% 
%     toc

    %% w visualization
    figure(1) ;
    
    w = out.w ;
    N = 100 ;
    W = reshape(full(msubs(w, [x;k], XK)),N,N) ;
    cc = contour(xvec, kvec, W) ; hold on ;
    clabel(cc);
    contour(xvec, kvec, W,[1 1],'LineWidth',1.5,'Color',plotting_color)
    
    % %% g visualization
    % tvec = linspace(0,1,500) ;
    % gvec = dmsubs(g,t,tvec) ;
    
    % figure(2) ; hold on ;
    % plot(tvec,gvec) 
    % plot([0 1],[1 1],'r')

    [grid_x, grid_y] = meshgrid(linspace(X_range(1, 1), X_range(1, 2), 50), linspace(K_range(1, 1), K_range(1, 2), 50));
    flat = [grid_x(:), grid_y(:)];
    spot_vals = reshape(dmsubs(out.w, [x; k], flat')', 50, 50);
    alf_vals = reshape(dmsubs(out_a.w, [x; k], flat')', 50, 50);

    COS(:,:,1) = zeros(50); % red
    COS(:,:,2) = zeros(50); % green
    COS(:,:,3) = zeros(50); % blue

    surf_a = surf(grid_x, grid_y, alf_vals, COS, 'FaceAlpha',0.5);
    surf_c = surf(grid_x, grid_y, spot_vals, 'FaceAlpha',0);
    surf_c.EdgeColor='red';
end
