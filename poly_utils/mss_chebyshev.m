function polys = mss_chebyshev(n, degree)
%degree is the maximum degree (inclusive)
%returns:
% polys: polys(n, l) is an mss poly representation of T_l(n)
zers = fliplr(eye(degree+1));
x = msspoly('x', n);
polys(1:n, 1:degree+1) = x(1);
for i=0:degree
    for j = 1:n
        polys(j, i+1) = flipud(cheb2mon(zers(:, i+1)))' * monomials(x(j), 0:degree);
    end
end

end