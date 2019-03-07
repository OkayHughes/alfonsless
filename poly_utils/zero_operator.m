function zero_arr = zero_operator(in_basis, out_basis)
    zero_arr = zeros(size(out_basis.monomials, 1), size(in_basis.monomials, 1));
end