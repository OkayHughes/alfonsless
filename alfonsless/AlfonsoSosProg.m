
% TODO:
%  -- test basis for dep. on indeterminates only.
%
classdef AlfonsoSosProg < handle
    properties (Access = protected)
        indeterminates = msspoly([]);
        dec_bases = [];

        constraint_polys = msspoly([]);

        dec_polys = msspoly([]);

        
    end

    properties (Abstract)
        K_wts;
    end
    

    
    methods (Access = protected)
        function res = exists_mon_basis(prog, vars, d)
            res=0;
            for mon_base_ind=1:size(prog.dec_bases, 1)
                mon_base = prog.dec_bases(mon_base_ind, 1);
                if var_set_eq(mon_base.variables, vars) && mon_base.d == d
                    res = mon_base_ind;
                end
            end
        end 

        function [As, bs, cs, const_bases] = monomial_formulate(prog, cost)
            bs = -1 * cost_vector_from_poly(cost, prog.dec_bases);
            As = {};
            cs = {};
            const_bases = [];
            for i=1:size(prog.constraint_polys, 1)
                [As_i, c_i, out_basis_i] = constraint_matrix_from_poly(prog.constraint_polys(i), ...
                                                      prog.dec_bases);
                As = [As; As_i];
                cs = [cs; c_i];
                const_bases = [const_bases; out_basis_i];
            end

        end

    end
    
    methods
        function prog = AlfonsoSosProg()

        end

        function [f, fcoeff] = new_free_poly(prog, vars, d)
            if ~var_set_subset(vars, prog.indeterminates))
                error('Variables given are not all registered with the alfonso program');
            end

            num_vars = size(vars, 1);
            sz = nchoosek(num_vars + d, d);
            mon_basis_ind  = prog.exists_mon_basis(vars, d);
            if mon_basis_ind
                mon_basis = prog.mononial_bases(mon_basis_ind); 
            else
                mon_basis = MonomialBasis(vars, d);
            end

            f_obj = Polynomial(mon_basis);

            prog.dec_bases = [prog.dec_bases; f_obj.mon_basis];

            f = transpose(f_obj.coeff_vars) * f_obj.mon_basis.monomials;
            fcoeff = f_obj.mon_basis.coeff_vars;
        end

        function with_indeterminate(prog, variables)
            prog.indeterminates = var_set_union(prog.indeterminates, variables);
        end
    

        
    end

    methods (Abstract)
        minimize(prog, cost);
        sos_on_K(prog, constraint_poly);
    end

end
