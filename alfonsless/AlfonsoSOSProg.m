
% TODO:
%  -- test basis for dep. on indeterminates only.
%
classdef AlfonsoSOSProg < handle
    properties (Access = protected)
        indeterminates = msspoly([]);
        dec_bases = [];

        constraint_polys = msspoly([]);
        constraint_degrees = [];
        constraint_variables = {};

        dec_polys Polynomial;
        
        
        A double;
        b double;
        c double;

        
    end


    

    
    methods (Access = protected)
        function res = exists_mon_basis(prog, vars, d)
            res=0;
            for mon_base_ind=1:size(prog.dec_bases, 1)
                mon_base = prog.dec_bases(mon_base_ind, 1);
                if var_set_eq(mon_base.variables, vars) && (mon_base.d == d)
                    res = mon_base_ind;
                end
            end
        end 

        function [As, bs, cs, const_bases] = monomial_formulate(prog, cost)
            bs = cost_vector_from_poly(cost, prog.dec_polys);
            As = {};
            cs = {};
            const_bases = [];
            for i=1:size(prog.constraint_polys, 1)
                [As_i, c_i, out_basis_i] = constraint_matrix_from_poly(prog.constraint_polys(i), ...
                                                      prog.dec_polys, prog.constraint_variables{i}, prog.constraint_degrees(i));
                As = [As; As_i];
                cs = [cs; c_i];
                const_bases = [const_bases; out_basis_i];
            end

        end

    end
    
    methods
        function prog = AlfonsoSOSProg()

        end

        function [f, fcoeff, f_monomial] = new_free_poly(prog, vars, d)
            if ~var_set_subset(vars, prog.indeterminates)
                error('Variables given are not all registered with the alfonso program');
            end

            mon_basis_ind  = prog.exists_mon_basis(vars, d);
            if mon_basis_ind
                mon_basis = prog.dec_bases(mon_basis_ind); 
            else
                mon_basis = MonomialBasis(vars, d);
            end

            assert(var_set_eq(vars, mon_basis.variables));

            f_obj = Polynomial(mon_basis);

            prog.dec_bases = [prog.dec_bases; f_obj.mon_basis];
            prog.dec_polys = [prog.dec_polys; f_obj];

            f = f_obj.to_msspoly();
            fcoeff = f_obj.coeff_vars;
            f_monomial = f_obj.mon_basis.monomials;
        end

        function with_indeterminate(prog, variables)
            prog.indeterminates = var_set_union(prog.indeterminates, variables);
        end
        
        
        function res = solve_alfonso(prog, g_h_params, in_opts)
            prob_data.A = prog.A;
            prob_data.b = prog.b;
            prob_data.c = prog.c;
            % make initial primal iterate
            g_h_params.U_arr;
            x0 = ones(sum(g_h_params.U_arr), 1);
            [~, g0, ~, ~] = alfonso_grad_and_hess(x0, g_h_params);
            
            rP = max((1+abs(prob_data.b))./(1+abs(prob_data.A*x0)));
            rD = max((1+abs(g0))./(1+abs(prob_data.c)));
            x0 = repmat(sqrt(rP*rD),sum(g_h_params.U_arr),1);
            if isfield(in_opts, "save_dir")
                opts.save_dir = fullfile(pwd(), in_opts.save_dir);
            end

            % run alfonso
            opts.optimTol = 1e-6 ;
            if isfield(in_opts, "verbose")
                opts.verbose = in_opts.verbose;
            end
            %opts.maxItRefineSteps = 2;
            
            res = alfonso(prob_data, x0, @alfonso_grad_and_hess_verbose, g_h_params, opts);
        end

        
    end

    methods (Abstract)
        minimize(prog, cost);
        sos_on_K(prog, constraint_poly);

    end

    properties (Abstract, Access = protected)
        K_wts;
    end

end
