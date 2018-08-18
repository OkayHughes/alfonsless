classdef MonomialBasis < handle
properties 
    n double;
    d double;
    num_monomials double;
    variables msspoly;
    monomials msspoly;
    power_matrix double;

end

methods
    function basis = MonomialBasis(variables, def_mons)

        if isa(variables, 'double')
            if ~isa(def_mons, 'double')
                error('second argument to MonomialBasis() is not numeric');
            end
            basis.variables = msspoly('@@@@', variables);
            basis.n = variables;
            basis.d = def_mons;
            basis.monomials = monomials(basis.variables, 0:basis.d);
        elseif isa(variables, 'msspoly')
            basis.variables = variables;
            basis.n = size(variables, 1);
            if isa(def_mons, 'double')
                basis.d = def_mons;
                basis.monomials = monomials(basis.variables, 0:basis.d);
            elseif isa(def_mons, 'msspoly')
                basis.monomials = def_mons;
                basis.d = max(arrayfun(msspoly_degree, def_mons));
            end
        end

        basis.num_monomials = size(basis.monomials, 1);

        basis.power_matrix = monomials_to_matrix(basis.monomials, basis.variables);
    end


end

end