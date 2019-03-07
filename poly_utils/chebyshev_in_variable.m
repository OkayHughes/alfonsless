function poly = chebyshev_in_variable(cheb_coeffs, dimension, basis)
    mon_coeffs = flipud(cheb2mon(cheb_coeffs));
    poly = mon_coeffs' * monomials(basis.variables(dimension), 0:(size(mon_coeffs, 1)-1));
end