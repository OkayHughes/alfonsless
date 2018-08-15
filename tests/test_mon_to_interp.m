degree = 6 ;
intParams = FeketeCube(2,degree) ;
basis = intParams.mon_basis;


% polynomial variable
t = basis.variables;

warn_flag = false;
err = 0;
for i=1:100
    vec = rand(size(basis.monomials));
    P0 = intParams.mon_to_P0 * vec;
    mon = vec' * basis.monomials;
    vals1 = dmsubs(mon, t, intParams.pts')';
    mon2 = (intParams.P0_to_mon * vals1)' * basis.monomials;
    vals2 = dmsubs(mon2, t, intParams.pts')';
    P0_to_vals2 = norm(P0 - vals2);
    vals1_to_vals2 = norm(vals1 - vals2);
    err = err +  max(P0_to_vals2, vals1_to_vals2);
    if (P0_to_vals2 > 1E-8 )| (vals1_to_vals2 > 1E-8)
        warn_flag = true;
    end
end


% for i=1:10
%     if mod(i, 2) == 1
%         mon_vec = rand(size(basis.monomials));
%         mon = mon_vec' * basis.monomials;
%         interp_vec = dmsubs(mon, basis.variables, intParams.pts')';
%     else
%         mon_vec = rand(size(lbasis.monomials));
%         mon = mon_vec' * lbasis.monomials;
%         interp_vec = dmsubs(mon, lbasis.variables, intParams.pts')';
%     end

%     mon = (intParams.P0_to_mon * interp_vec)' * basis.monomials;
%     evals = dmsubs(mon, basis.variables, intParams.pts')';
%     normm = norm(evals-interp_vec);
%     if normm > 1E-10 & mod(i, 2) == 1
%         warn_flag = true;
%     end
% end

if warn_flag
    warning('P0_to_mon may be badly behaved\n err = %d, should be very small', err)
end