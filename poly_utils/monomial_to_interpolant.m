function [P0_to_mon, mon_to_P0] = monomial_to_interpolant(interpolant_basis, mss_chebyshev, basis)
    %Takes:
    %   interpolant_basis: U X L real matrix
    %       interpolant basis as calculated in FeketeCube
    %   interpolant_basis: U X L real matrix
    %       Q in QR decomposition of interpolant_basis
    %   mss_chebyshev: L x 1 msspoly matrix
    %       mss_chebyshev(i) is the polynomial corresponding to
    %       interpolant_basis(i)
    %   basis: struct
    %       Struct returned from `monomial_basis`
    %Returns:
    %   mon_to_interp: U x L real matrix
    %       converts polynomial in monomial basis (specified in basis) 
    %       to interpolant basis
    %   interp_to_mon: L x U real matrix 
    %       converts polynomial in interpolant basis to the monomial basis
    %       specified in basis
    
    M = zeros(size(mss_chebyshev, 1));
    for mon_ind=1:size(mss_chebyshev, 1)
        M(:, mon_ind) = msspoly_to_vector(mss_chebyshev(mon_ind), basis);
    end
    cond_num = cond(full(M));
    M_inv = inv(M);
    if cond_num > 1E10
        warning('Condition number of polynomial basis is large (cond(M) = %d)\nCOB matrices may be inaccurate.', cond_num);
    end
    cond_num = cond(full(interpolant_basis));
    if cond_num > 1E10
        warning('Condition number of P0 is large (cond(P0) = %d)\nCOB matrices may be inaccurate.', cond_num);
    end
    
    P0_to_mon = M * inv(interpolant_basis);
    mon_to_P0 = interpolant_basis * M_inv;
end

