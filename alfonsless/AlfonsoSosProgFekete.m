
% TODO:
%  -- test basis for dep. on indeterminates only.
%
classdef AlfonsoSosProgFekete < handle
    properties (Access = protected)

        K_wts = {};

    end
    
    methods
        function prog = AlfonsoSosProgFekete()

        end
        
        function sol = minimize(prog, cost, bounds)
            [As, bs, cs, const_bases] = prog.monomial_formulate(cost);
            dec_interp_bases = [];
            dec_mon_to_interps = {};
            dec_interp_to_mons = {};

            for i=1:size(prog.dec_polys, 1)
                pol = prog.dec_polys(i);

                dec_interp_basis = AffineFeketeCube(pol.mon_basis.variables, ...
                                                    pol.mon_basis.d, bounds{i});
                dec_interp_bases = [dec_interp_bases; dec_interp_basis];

                [dec_interp_to_mon, dec_mon_to_interp] = dec_interp_basis.inter_to_mon(prog.dec_polys(i).mon_basis);
                dec_mon_to_interps = [dec_mon_to_interps; dec_mon_to_interp];
                dec_interp_to_mons = [dec_interp_to_mons; dec_interp_to_mon];

            end
            b = [];

            for i=1:size(bs, 1)
                b_i = bs(i)' * dec_interp_to_mons{i};
                b = [b; -b_i]; 
            end

            const_interp_bases = [];
            const_mon_to_interps = {};
            const_interp_to_mons = {};
            for i=1:size(const_bases, 1)
                const_interp_basis = AffineFeketeCube(const_bases(i).variables, ...
                                                      const_bases(i).d, ...
                                                      prog.K_wts{i});
                const_interp_bases = [const_interp_bases; const_interp_basis];
                [const_interp_to_mon, const_mon_to_interp] = const_interp_basis.inter_to_mon(const_bases(i));
                const_mon_to_interps = [const_mon_to_interps; const_mon_to_interp];
                const_interp_to_mons = [const_interp_to_mons; const_interp_to_mon];
            end


            for i=1:size(As, 1)
                for j=1:size(As, 2)
                    As{i} = -const_mon_to_interps{i} * As{i, j} * const_interp_to_mons{j};
                    cs{i} = -const_mon_to_interps{i} * cs{i};
                end
            end
            A = cell2mat(As);
            c = cell2mat(cs);

            interp_params_arr = const_interp_bases';
            g_h_params = gen_grad_params(interp_params_arr);

            % create A,b,c matrices to define conic problem
            prob_data.A = A;
            prob_data.b = b;
            prob_data.c = c;

            % make initial primal iterate

            x0 = ones(sum(g_h_params.U_arr), 1);
            [~, g0, ~, ~] = alfonso_grad_and_hess(x0, g_h_params);

            rP = max((1+abs(prob_data.b))./(1+abs(prob_data.A*x0)));
            rD = max((1+abs(g0))./(1+abs(prob_data.c)));
            x0 = repmat(sqrt(rP*rD),sum(g_h_params.U_arr),1);


            % run alfonso
            opts.optimTol = 1e-6 ;
            results = alfonso(prob_data, x0, @alfonso_grad_and_hess, g_h_params, opts);

            res_polys = msspoly(zeros(size(prog.dec_vars, 1)));
            sm = 1;
            for i=1:size(prog.dec_vars, 1)
                res_polys(i) = dec_interp_to_mons * results.y(sm:sm - 1 + g_h_params.U_arr(i), 1);
                sm = sm + g_h_params.U_arr(i);
            end

            sol = res_polys;


        end
    
        function sos_on_K(prog, constraint, bounds)
            prog.constraint_polys = [prog.constraint_polys; constraint];
            prog.K_wts = [prog.K_wts; bounds];
        end
        
    end


end
