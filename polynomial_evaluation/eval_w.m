function w_out = eval_w(k,wcoeffs,wpowers)
% w_out = eval_w(k,wcoeffs,wpowers)
%
% Evaluate a polynomial w(k), expressed as a coefficient matrix and a
% powers matrix, at the point k (2-by-1)
    w_out = wcoeffs*prod(repmat(k,1,size(wpowers,1))'.^wpowers,2) ;
end