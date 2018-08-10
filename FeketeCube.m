% This code generates parameters for the interpolant basis representation
% of sum-of-squares polynomials with an arbitrary number of variables.
% It follows the approach described in:
%
% A. Sommariva and M. Vianello, Computing approximate Fekete points by QR
% factorizations of Vandermonde matrices, Computers & Mathematics with
% Applications, 57 (2009), pp. 1324-1336. Available at
% https://doi.org/10.1016/j.camwa.2008.11.011.
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Date: 06/14/2018
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
%
% partitions. The partitions function computes all partitions of an integer.
% We use the implementation of John D'Errico from a MATLAB Central File 
% Exchange post which is available at
% https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer.
% -------------------------------------------------------------------------


function intParams = FeketeCube(n,d, variables)
% This method generates parameters for the interpolant basis representation
% of sum-of-squares polynomials with an arbitrary number of variables.
% --------------------------------------------------------------------------
% USAGE of "FeketeCube"
% intParams = FeketeCube(n,d)
% --------------------------------------------------------------------------
% INPUT
% n:	number of arguments to the polynomials
% d:	degree of polynomials to be squared
%
% OUTPUT
% intParams:                interpolation parameters
% - intParams.n:            number of arguments to the polynomials
% - intParams.d:            degree of polynomials to be squared
% - intParams.L:            dimension of the space of (intParams.n)-variate
%                           degree-d polynomials
% - intParams.U:            dimension of the space of (intParams.n)-variate
%                           degree-(2*d) polynomials
% - intParams.pts:          approximate Fekete points for degree-(2*d)
%                           polynomial interpolation. (intParams.U x 1) array.
% - intParams.w:            (scaled) weights for Clenshaw-Curtis quadrature
% - intParams.P0:           evaluations of bivariate product Chebyshev 
%                           polynomials of the first kind up to degree d at
%                           the points intParams.pts
%                           (intParams.U x intParams.L) array.
% - intParams.P:            evaluations of a basis for the space of
%                           (intParams.n)-variate degree-d polynomials 
%                           at the points intParams.pts. 
%                           (intParams.U x intParams.L) array with 
%                           orthonormal columns.
% - intParams.nrPoints1D:	smallest number of points along a side of the
%                           initial interpolation grid
% - intParams.nrPoints:     number of points in the initial interpolation 
%                           grid
% - intParams.mon_basis     basis object returned by monomial_basis
% - intParams.P0_to_mon change-of-basis matrix taking vectors expressed in
%                           the columns of P0 into the monomial basis specified
%                           by mon_basis
% - intParams.P_to_mon change-of-basis matrix taking vectors expressed in
%                           the columns of P into the monomial basis specified
%                           by mon_basis
% - intParams.mon_to_P0     inverse of P0_to_mon(in some sense)
% - intParams.mon_to_P      inverse of P_to_mon(in some sense)    
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% chebpolyval from Chebfun. Chebfun is an open-source package for  
% numerical computing with functions: http://www.chebfun.org/.
% Our code has been developed and tested with Chebfun version 5.5.0.
% The latest version of the Chebfun package can be downloaded from
% http://www.chebfun.org/download/.
%
% partitions. The partitions function computes all partitions of an integer.
% We use the implementation of John D'Errico from a MATLAB Central File 
% Exchange post which is available at
% https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer.
% -------------------------------------------------------------------------
  
    intParams.n = n;
    intParams.d = d;    
    intParams.L = nchoosek(n+d,n);
    intParams.U = nchoosek(n+2*d,n);
    if nargin < 3
        intParams.mon_basis = monomial_basis(intParams.n, intParams.d*2, 'f');
    else
        intParams.mon_basis = monomial_basis(intParams.n, intParams.d*2, variables);
    end
    
    intParams.nrPoints1D = 2*d+1;

    intParams.nrPoints = prod(intParams.nrPoints1D:intParams.nrPoints1D+n-1);
    pts = zeros(intParams.nrPoints,n);
    for j = 1:n
        temp = 1;
        for i = 1:j-1; temp = kron(temp,ones(intParams.nrPoints1D+i-1,1)); end;
        temp = kron(temp,chebpts(intParams.nrPoints1D+j-1));
        for i = j+1:n; temp = kron(temp,ones(intParams.nrPoints1D+i-1,1)); end;
        pts(:,j) = temp;
    end
    
    P = ones(intParams.nrPoints,intParams.U);
    m = ones(intParams.U,1);
    
    prod_polynomials = msspoly(ones(intParams.U, 1));
    col     = 0;
    lrEye   = fliplr(eye(2*d+1)); %descending T_n
    for t = 0:2*d      % polynomials with total degree up to 2*d
        allDegs = partitions(t, ones(1,n));
        [nrDegs,~] = size(allDegs);
        for i = 1:nrDegs
            col = col+1;
            for j = 1:n
                dj = allDegs(i,j);
                prod_polynomials(col) = prod_polynomials(col) * ...
                                            chebyshev_in_variable(lrEye(:,dj+1), ...
                                                                  j, ...
                                                                  intParams.mon_basis);
                P(:,col) = P(:,col).*chebpolyval(lrEye(:,dj+1),pts(:,j));


                if dj == 1; m(col) = 0;
                else; m(col) = m(col)*(((-1)^dj+1)/(1-dj^2)); end;
            end
        end
    end
    
    w = P'\m;
    ind = abs(w)>0;
    
    % extracts the positive entries of w
    w = w(ind);
    % extracts the subset of points indexed with the support of w
    pts = pts(ind,:);
    % extracts the subset of polynomials up to total degree d
    P0_large = P(ind, :);
    P = P(ind,1:intParams.L);


    
    %P' should equal P_test
    P_test = dmsubs(prod_polynomials, intParams.mon_basis.variables, pts');
    normm = norm(P_test - P0_large', 'fro');
    if norm(P_test - P0_large', 'fro') > 1E-8
        warning('Interpolant basis monomials may not be accurate\n ||P0 - P_test|| = %d, should be very small', normm);
    end
    intParams.w = w;
    intParams.pts = pts;
    intParams.P0 = P;
    
    [P_large, ~] = qr(P0_large);
    [intParams.P,~] = qr(P,0);
    
    intParams.P_full = P_large;
    intParams.P0_full = P0_large;
    
    [P0_to_mon, P_to_mon, mon_to_P0, mon_to_P] = monomial_to_interpolant(P0_large, P_large, prod_polynomials, intParams.mon_basis);
    
    intParams.polynomials = prod_polynomials;
    
    intParams.P0_to_mon = P0_to_mon;
    intParams.P_to_mon = P_to_mon;
    intParams.mon_to_P0 = mon_to_P0;
    intParams.mon_to_P = mon_to_P;
    
      
return
    