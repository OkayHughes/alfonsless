%% user parameters
degree = 8 ;
solver = 'mosek' ;
visualize = 1 ;
initial_cond_box_side_width = 0.2 ;

k_viz = [0.5 0 -0.5] ; % for visualizing w > 1 at different values of k;
                       % can be vector or scalar

t_viz = 0:0.25:1 ; % for visualizing the "time slices" in the dec variable
                   % v, which captures the flow of the dynamics f + g; this
                   % will be done for the last element of k_viz

%% set up advection vars
t = msspoly('t', 1);
z = msspoly('z', 2);
k = msspoly('k', 1) ;

%% set up spaces
T = 1 ;

Z_range = repmat([-1, 1], size(z)) ;

L = initial_cond_box_side_width/2 ; 
Z0_range = repmat([-L, L], size(z)) ;

K_range = repmat([-1,1], size(k)) ;

%% dynamics
% system dynamics
f = [0.5*k ;
     0.25*k ] ;

% error dynamics
g = [0 ;
     0.25*k*t] ;

%% create domain
hZ = (z - Z_range(:,1)).*(Z_range(:,2) - z) ;

hZ0 = (z - Z0_range(:,1)).*(Z0_range(:,2) - z) ;
% hX0 = 1 - (x(1)/X0_range(1,1))^2 - (x(2)/X0_range(1,2))^2 ;

hK = (k - K_range(:,1)).*(K_range(:,2) - k) ;

%% create cost function
int_ZK = boxMoments([z;k], [Z_range(:,1); K_range(:,1)], [Z_range(:,2); K_range(:,2)]) ;

%% setup the problem structure
prob.t = t ;
prob.z = z ;
prob.k = k ;
prob.f = f ;
prob.g = g ;
prob.hZ = hZ ;
prob.hZ0 = hZ0 ;
prob.hK = hK ;
prob.int_ZK = int_ZK ;
prob.degree = degree ;
prob.domain_size = 1 ;
prob.solver = solver ;
prob.run_solver = 1 ;

%% solve the problem
out = findFRS(prob) ;

prob_a.degree = degree;
prob_a.t = t;
prob_a.x = z;
prob_a.k = k;
prob_a.f = f;
prob_a.g = g;
prob_a.X_bounds = Z_range;
prob_a.X0_bounds = Z0_range;
prob_a.K_bounds = K_range;
prob_a.T = T;

out_a = find_FRS_alfonso(prob_a);

Y_bounds = [Z_range; K_range];
l2_dist = l2_dist_on_box(out.w, out_a.w, Y_bounds(:, 1), Y_bounds(:, 2), [z;k])

%% set up variables for visualization
N = 100 ;

xvec = linspace(Z_range(1,1),Z_range(1,2),N) ;
yvec = linspace(Z_range(2,1),Z_range(2,2),N) ;

XY = makeContourAxes(N) ;
[X,Y] = makeContourAxes(N) ;


for kval = k_viz
    figure('Name',sprintf('w_a vs w_s k=%d', kval)) ; cla ; hold on ;
    [grid_x, grid_y] = meshgrid(linspace(Z_range(1, 1), Z_range(1, 2), 50), linspace(Z_range(2, 1), Z_range(2, 2), 50));
    flat = [grid_x(:), grid_y(:), kval * ones(50^2, 1)];
    spot_vals = reshape(dmsubs(out.w, [z; k], flat')', 50, 50);
    alf_vals = reshape(dmsubs(out_a.w, [z; k], flat')', 50, 50);

    COS(:,:,1) = zeros(50); % red
    COS(:,:,2) = zeros(50); % green
    COS(:,:,3) = zeros(50); % blue

    surf_a = surf(grid_x, grid_y, alf_vals, COS, 'FaceAlpha',0.5);
    surf_c = surf(grid_x, grid_y, spot_vals, 'FaceAlpha',0);
    surf_c.EdgeColor="red";
end


%% output visualization
if visualize
    figure(1) ; clf ;
    
    %% w visualization
    subplot(1,2,1) ; hold on ; axis equal
    
    w = out.w ;
    
    for kidx = k_viz
        wk = msubs(w, k, kidx) ;
        W = reshape(full(msubs(wk,z,XY)),N,N) ;
        contour(X,Y,W,[1 1],'LineWidth',1.5,'Color',[0 0.8 0.2])
    end
    
    X0 = reshape(full(msubs(hZ0(1),z(1),XY(1,:))),N,N) ;
    Y0 = reshape(full(msubs(hZ0(2),z(2),XY(2,:))),N,N) ;
    Z0 = (X0 >= 0) & (Y0 >= 0) ;
    contour(X,Y,Z0,[1 1],'LineWidth',1.5,'Color',[0 0 1])
    
    title('w >= 1 for variety of k')
    
    % %% v visualization
    % subplot(1,2,2) ; hold on ; axis equal
    % disp('Visualizing v for the last k in k_for_viz')
    
    % v = out.v ;
    % kidx = k_viz(end) ;
    % vk = msubs(v,k,kidx) ;
    
    % for tidx = t_viz
    %     vt = msubs(vk,t,tidx) ;
    %     V = reshape(full(msubs(vt,z,XY)),N,N) ;
    %     contour(X,Y,V,[0 0],'LineWidth',1.5,'Color',[0 0.8 0.2])
    % end
    
    % X0 = reshape(full(msubs(hZ0(1),z(1),XY(1,:))),N,N) ;
    % Y0 = reshape(full(msubs(hZ0(2),z(2),XY(2,:))),N,N) ;
    % Z0 = (X0 >= 0) & (Y0 >= 0) ;
    % contour(X,Y,Z0,[1 1],'LineWidth',1.5,'Color',[0 0 1])
    
    % title('v <= 0 at varying time t')
end
