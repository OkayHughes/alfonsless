function J_out = eval_J(k,Jcoeffs,Jpowers,N)
% J_out = eval_J(k,Jcoeffs,Jpowers,N)
%
% Given a Jacobian matrix represented as a coefficients matrix and a powers
% matrix (see diff_wk_wrt_k), evaluate the jacobian at the point k (2-by-1)

    J_out = zeros(N,2) ;
    J_out(:,1) = Jcoeffs(:,1:N)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,1:2),2)) ;
    J_out(:,2) = Jcoeffs(:,N+1:end)'*(prod(repmat(k',size(Jpowers,1),1).^Jpowers(:,3:4),2)) ;
end