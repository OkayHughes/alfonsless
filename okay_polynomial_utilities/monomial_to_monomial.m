function [arg1_to_arg2, arg2_to_arg1] = monomial_to_monomial(basis1, basis2)
    %assumes that one basis has greater variables, and dominates the degree of the other
    %bases also must be written in the same msspoly free variables
    %operates on the general principle that 
    %(arg1_to_arg2 * msspoly_to_vector(f, basis1))' * basis2.monomials = f

    %flip order so I can just iterate over the longer set of monomials
    [~, max_ind] = max([size(basis1.monomials, 1), size(basis2.monomials, 1)]);
    if max_ind == 2
        long = zeros(size(basis2.monomials, 1), size(basis1.monomials, 1));
        short = zeros(size(basis1.monomials, 1), size(basis2.monomials, 1));

        sbasis = basis1;
        lbasis = basis2;
    else
        short = zeros(size(basis2.monomials, 1), size(basis1.monomials, 1));
        long = zeros(size(basis1.monomials, 1), size(basis2.monomials, 1));

        sbasis = basis2;
        lbasis = basis1;
    end

    for i=1:size(lbasis.monomials, 1)
        if i <= size(sbasis.monomials, 1)
            vec_l = msspoly_to_vector(sbasis.monomials(i), lbasis);
            long(:, i) = vec_l;
        end
        try
            vec_s = msspoly_to_vector(lbasis.monomials(i), sbasis);
            short(:, i) = vec_s;
        catch
        end
    end



    if max_ind == 2
        arg1_to_arg2 = long;
        arg2_to_arg1 = short;
    else
        arg1_to_arg2 = short;
        arg2_to_arg1 = long;
    end
end