function zd = segway_dynamics(t,z,k,hdmax,vmax,braking)
    if nargin < 6
        braking = false ;
    end
    
    % states
    h = z(3) ;
    hd = z(4) ;
    v = z(5) ;
    
    % inputs
    hddes = hdmax*k(1) ;
    Kg = 2.95 ;
    g = Kg*(hddes - hd) ;
    
    vdes = (vmax/2)*k(2) + vmax/2 ;
%     vdelta = vdes - v ;
%     if vdelta > 0
%         Ka = 3 ;
%     else
%         Ka = 9.81/2 ;
%     end
    Ka = 3 ;
    
    if braking && t > 0.5
        a = Ka*(-v) ;
    else
        a = Ka*(vdes - v) ;
    end
    
    zd = [v*cos(h) ;
          v*sin(h) ;
          hd ;
          g ;
          a ] ;
end