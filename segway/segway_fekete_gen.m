function segway_FRS_solver_okay(in_degree)
setup()
% options
degree = in_degree ;
solver = 'mosek' ;
notes = '' ;
run_solver = 1 ;

filename = ['segwayFRS_deg', num2str(degree), notes,'.mat'] ; 

% initial conditions
x0 = -0.75 ;
y0 = 0 ;

% vehicle parameters
rcar = 0.2 ;
vmax = 0.5 ; % vmax = 2 usually
wmax = 1 ;

% space parameters
T = 1 ;
D = 1.5 ;

%% set up advection vars
t = msspoly('t', 1);
z = msspoly('z', 2);
k = msspoly('k', 2) ;

x = z(1) ;
y = z(2) ;

kw = k(1) ;
kv = k(2) ;

w = kw*wmax ;
v = (vmax/2)*kv + (vmax/2) ;

%% set up spaces
Z_range = repmat([-1, 1], size(z)) ;
Z0_range = 1/D * [x0 - rcar, x0 + rcar;
                  y0 - rcar, y0 + rcar ];

K_range = repmat([-1,1], size(k)) ;

%% dynamics
% system dynamics
f = (T/D)*[v - w*(D*y-y0) ;
           w*(D*x-x0) ] ;

g = (T/D)*[(t-1).^2 , 0 ;
           0 , (D*x-x0)*(t-1).^2 ] ;

%% create domain
hZ = (z - Z_range(:,1)).*(Z_range(:,2) - z) ;

%circ_hZ0 = 1 - ((x-x0/D)/(rcar/D)).^2 - ((y-y0/D)/(rcar/D)).^2 ;

hZ0 = (z - Z0_range(:,1)).*(Z0_range(:,2) - z) ;

hK = (k - K_range(:,1)).*(K_range(:,2) - k) ;


%% create cost function
int_ZK = boxMoments([z;k], [Z_range(1:2,1); K_range(:,1)], [Z_range(1:2,2); K_range(:,2)]) ;

%% setup the problem structure
% prob.t = t ;
% prob.z = z ;
% prob.k = k ;
% prob.f = f ;
% if exist('g','var')
%     prob.g = g ;
% end
% prob.hZ = hZ ;
% prob.hZ0 = hZ0 ;
% prob.hK = hK ;
% prob.T = 1 ;
% prob.int_ZK = int_ZK ;
% prob.degree = degree ;
% prob.domain_size = 1 ;
% prob.solver = solver ;
% prob.run_solver = run_solver ;


% %% solve the problem
% out = findFRS(prob) ;

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
prob_a.verbose = true;

out_a = gen_fekete_FRS(prob_a);

% Y_bounds = [Z_range; K_range];
% l2_error_per_unit = l2_dist_on_box(out.w, out_a.w, Y_bounds(:, 1), Y_bounds(:, 2), [z;k]) / prod(Y_bounds(:, 2) - Y_bounds(:, 1))



%save(filename)
end
