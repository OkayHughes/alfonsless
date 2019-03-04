
classdef AlfonsoSOSProgFekete < AlfonsoSOSProg & handle
    properties (Access = protected)

        K_wts = {};

    end
    
    methods
        function prog = AlfonsoSOSProgFekete()

        end
        
        function sol = minimize(prog, cost)
            bounds = cell(size(prog.dec_polys));
            for i=1:size(prog.dec_polys, 1)
                bounds{i} = repmat([-1, 1], prog.dec_polys(i).mon_basis.n, 1);
            end
            
            [As, bs, cs, const_bases] = prog.monomial_formulate(cost);
            dec_interp_bases = [];
            dec_mon_to_interps = {};
            dec_interp_to_mons = {};

            for i=1:size(prog.dec_polys, 1)
                pol = prog.dec_polys(i);
                dec_interp_basis = AffineFeketeBasis(pol.mon_basis.variables, ...
                                                    pol.mon_basis.d, bounds{i});
                dec_interp_bases = [dec_interp_bases; dec_interp_basis];

                [dec_interp_to_mon, dec_mon_to_interp] = dec_interp_basis.inter_to_mon(prog.dec_polys(i).mon_basis);
                dec_mon_to_interps = [dec_mon_to_interps; dec_mon_to_interp];
                dec_interp_to_mons = [dec_interp_to_mons; dec_interp_to_mon];

            end
            b = [];

            for i=1:size(bs, 1)
                b_i = bs{i}' * dec_interp_to_mons{i};
                b = [b; -b_i']; 
                size(b_i)
            end
            

            const_interp_bases = [];
            const_mon_to_interps = {};
            const_interp_to_mons = {};
            for i=1:size(const_bases, 1)
                const_interp_basis = AffineFeketeBasis(const_bases(i).variables, ...
                                                      const_bases(i).d, ...
                                                      prog.K_wts{i});
                const_interp_bases = [const_interp_bases; const_interp_basis];
                [const_interp_to_mon, const_mon_to_interp] = const_interp_basis.inter_to_mon(const_bases(i));
                const_mon_to_interps = [const_mon_to_interps; const_mon_to_interp];
                const_interp_to_mons = [const_interp_to_mons; const_interp_to_mon];
            end
            for i=1:size(As, 1)
                for j=1:size(As, 2)
                    As{i, j} = -const_mon_to_interps{i} * As{i, j} * dec_interp_to_mons{j};
                end
                cs{i} = const_mon_to_interps{i} * cs{i};
            end
            
            A = cell2mat(As);
            c = cell2mat(cs);

            interp_params_arr = const_interp_bases';
            g_h_params = gen_grad_params(interp_params_arr);

            size(A);
            size(b);
            size(c);
            % create A,b,c matrices to define conic problem
            prog.A = A';
            prog.b = b;
            prog.c = c;
            results = prog.solve_alfonso(g_h_params);

            res_polys = msspoly(zeros(size(prog.dec_polys, 1), 1));
            sm = 1;
            for i=1:size(prog.dec_polys, 1)
                res_polys(i) = prog.dec_polys(i).eval(dec_interp_to_mons{i} * results.y(sm:sm - 1 + prog.dec_polys(i).mon_basis.num_monomials, 1));
                sm = sm + prog.dec_polys(i).mon_basis.num_monomials;
            end

            sol.polys = res_polys;
            sol.A = prog.A;
            sol.b = prog.b;
            sol.c = prog.c;
            
            sol.cost = -prog.b' * results.y;


        end
    
        function sos_on_K(prog, constraint, vars, bounds, degree)
            if mod(degree, 2) ~= 0
                error('Degree should be even')
            end
            if msspoly_degree(constraint) - 1 > degree
                error("Degree should exceed terms appearing in the constraint (requirement may go away in later versions)")
            end
            if size(bounds, 1) ~= size(vars, 1)
                error('Size of bounds should have the same number of rows as variables')
            end
            prog.constraint_polys = [prog.constraint_polys; constraint];
            prog.K_wts = [prog.K_wts; bounds];
            prog.constraint_degrees = [prog.constraint_degrees; degree];
            prog.constraint_variables = [prog.constraint_variables; {vars}];
        end
        
    end


end
