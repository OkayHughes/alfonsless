classdef Polynomial
properties 
    mon_basis MonomialBasis;
    coeff_vars msspoly;

end

methods
    function pol = MonomialBasis(mon_basis, coeff_vars)
        if nargin > 1
            pol.mon_basis = mon_basis;
            pol.coeff_vars = coeff_vars;
        else
            pol.mon_basis = mon_basis;
            pol.coeff_vars = Polynomial.new_coeff_vars(mon_basis);
        end
    end

    function pol = eval(in_pol, coeffs)
        pol = transpose(in_pol.coeff_vars) * in_pol.mon_basis.monomials;
        pol = msubs(pol, in_pol.coeff_vars, coeffs);
    end

end

methods (Static)
    function coeff_vars = new_coeff_vars(arg)
        id = MonomialBasis.new_unique_id();
        if isa(arg, 'msspoly')
            coeff_vars = msspoly(sprintf('@%s', id), size(arg, 1));
        elseif isa(arg, 'double')
            coeff_vars = msspoly(sprintf('@%s', id), arg);
        elseif isa(arg, 'MonomialBasis')
            coeff_vars = msspoly(sprintf('@%s', id), arg.num_monomials);
        end
    end
end


methods (Static, Access=protected)
    function un_id = new_unique_id()
        persistent COEFF_VAR_COUNTER;

        chars = 'abcdefghijklmnopqrstuvwxyz';

        if isempty(COEFF_VAR_COUNTER)
             COEFF_VAR_COUNTER = 0;
        end

        if COEFF_VAR_COUNTER < 26^3
            COEFF_VAR_COUNTER = COEFF_VAR_COUNTER + 1;
        else
            COEFF_VAR_COUNTER = 1;
            warning('coeff variables may not be unique');
        end
        un_id = '';
        tmp = COEFF_VAR_COUNTER;
        for pow=1:3
            divisor = mod(tmp, 26^pow);
            dividend = (26^(pow-1));
            amt1 = divisor/dividend;
            un_id = [un_id, chars(amt1+1)];
            tmp = tmp - divisor;
            if tmp==0
                break
            end

        end


    end
end

end