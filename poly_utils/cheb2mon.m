function b = cheb2mon(a)
%CHEB2MON  Chebyshev to monomial basis conversion.
%   B = CHEB2MON(A) converts polynomial A given in Chebyshev basis to 
%   monomial basis B. The polynomial must be given with its coefficients in
%   descending order, i.e. A = A_N*T_N(x) + ... + A_1*T_1(x) + A_0*T_0(x)
%
%   Example: 
%    Suppose we have a polynomial in Chebyshev basis: 
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x), where T_0=1, T_1=x, T_2=2x^2-1
%    and for example a2=1, a1=0, a0=-1.
%    We want to express the polynomial in the monomial base {1,x,x^2), i.e.
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x) = b2*x^2 + b1*x + b0,
%    where b = [b2 b1 b0] is sought.
%    Solution:
%      a = [1 0 -1];
%      b = cheb2mon(a);
%
%   See also   CHEB2MON, CHEBPOLY

%   Zoltán Csáti
%   2015/03/31


% Construct the Chebyshev polynomials of the first kind
N = length(a)-1;
C = chebpoly(N);
% Create the transformation matrix
A = zeros(N+1);
for k = 1:N+1
   A(k:N+1,k) = C{N+2-k};
end
% Perform the basis conversion using a matrix-vector product
b = A*a(:);


function C = chebpoly(n)
% C = CHEBPOLY(N) returns Chebyshev polynomials of the first kind from 
% degree 0 to n

C = cell(n+1,1);   % using cell array because of the different matrix sizes
C{1} = 1;          % T_0 = 1
C{2} = [1 0];      % T_1 = x
for k = 3:n+1      % T_n = 2xT_n-1 - T_n-2
    C{k} = [2*C{k-1} 0] - [0 0 C{k-2}];
end