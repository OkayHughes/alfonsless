function [X,Y] = makeContourAxes(N,B)
    if nargin < 1
        N = 100 ;
    elseif nargin < 2
        B = [-1,1] ;
    end
    
    v = linspace(B(1),B(2),N) ;
    
    [X,Y] = meshgrid(v,v) ;
    
    if nargout == 1
        X = [X(:), Y(:)]' ;
    end
end