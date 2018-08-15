%% setup
clear ; clc ;

load('segwayFRS_deg12_vmax1p5_D1p5_withUnc.mat')

w = out.w ;

k1vec = 0 ;
k2vec = -1 ;

NP = 0 ;

visualize_v_flag = false ;
dt = 0.25 ;
tvec = [0 0.8] ;

%% set up spaces
N = 100 ;

xvec = linspace(Z_range(1,1),Z_range(1,2),N) ;
yvec = linspace(Z_range(2,1),Z_range(2,2),N) ;

kxvec = linspace(K_range(1,1),K_range(1,2),N) ;
kyvec = linspace(K_range(2,1),K_range(2,2),N) ;

[X,Y] = meshgrid(xvec, yvec) ;
ZZ = [X(:) Y(:)]' ;

[KX,KY] = meshgrid(kxvec, kyvec) ;
KK = [KX(:) KY(:)]' ;


P = [0.1+0.4.*rand(1,NP) ;
     0 + 0.5.*(2*rand(1,NP)-1)] ;

colors = 0.7.*rand(NP,3) ;

%% w visualization
figure(223) ; clf ;
subplot(1,2,2)
hold on
tic

% plot segway
thetavec = linspace(0,2*pi,100) ;
segxvec = (rcar/D).*cos(thetavec);
segyvec = (rcar/D).*sin(thetavec);
plot(segxvec+x0/D,segyvec+y0/D) ;

for k1idx = k1vec
    wk1 = msubs(w,k(1),k1idx) ;
    for k2idx = k2vec
        wk2 = msubs(wk1,k(2),k2idx) ;

        F = reshape(full(msubs(wk2, [x;y], ZZ)),N,N) ;
        
        contour(xvec,yvec,F-1,[0 0],'Color',...
               [0, 0.7*abs(k1idx), 0.7*abs(k2idx)],'LineWidth',1.5)
        
        ktemp = [k1idx;k2idx] ;
           
        [tout,zout] = ode45(@(tau,zed) segway_dynamics(tau,zed,ktemp,D,hdmax,vmax),...
                            [0,1],[x0/D y0/D 0 0 vmax/D]') ;
            
        plot(zout(:,1),zout(:,2))
        
        if visualize_v_flag
            % v visualization
            for tidx = tvec
                vtemp = msubs(out.v,[t;k],[tidx;k1idx;k2idx]) ;
                F2 = reshape(full(msubs(vtemp,[x;y],ZZ)),N,N) ;
                contour(xvec,yvec,F2-1,[0,0], 'Color',...
                    [0, abs(k1idx), abs(k2idx)],'LineWidth',1.5)
            end
        end
        
        % plot car at end of trajectory
        plot(segxvec+zout(end,1),segyvec+zout(end,2),'b','LineWidth',1.5)
    end
end

axis equal
toc

%% K visualization
subplot(1,2,1)

tic
for pidx = 1:size(P,2)
%     random points
    xidx = P(:,pidx) ;
    
    widx = msubs(w,[x;y],xidx) ;
    Widx = reshape(full(msubs(widx,k,KK)),N,N) ;
    
    col = colors(pidx,:) ;
    
    subplot(1,2,2)
    hold on
    plot(xidx(1),xidx(2),'x','Color',col,'MarkerSize',15)
    axis equal
    
    subplot(1,2,1)
    hold on
    contour(kxvec,kyvec,Widx'-1,[0 0],'Color',col,'LineWidth',1.25)
    ylabel('steering')
    xlabel('speed')
    axis equal
end
toc