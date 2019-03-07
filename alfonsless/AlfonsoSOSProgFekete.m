
classdef AlfonsoSOSProgFekete < AlfonsoSOSProg & handle
    properties (Access = protected)

        K_wts = {};

    end
    
    methods
        function prog = AlfonsoSOSProgFekete()

        end
        
        function sol = minimize(prog, cost, opts)
            if isfield(opts, "verbose")
                verbose = opts.verbose;
            else
                verbose = 1;
            end
            
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
                                                    pol.mon_basis.d, bounds{i}, verbose);
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
                                                      prog.K_wts{i}, verbose);
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
            if exist("opts")
                results = prog.solve_alfonso(g_h_params, opts);
            else
                results = prog.solve_alfonso(g_h_params, struct());
            end

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
        
        %detail level 0 prints constraint degrees and dimensions
        %detail level 1 analyzes matrix properties in monomial space
        %detail level 2 analyzes matrix properties in interpolant space
        
        function problem_chars(prog, detail_level)
            
            if detail_level > 0
                [As, ~, ~, const_bases] = prog.monomial_formulate(msspoly(0));
                dec_interp_bases = [];
                dec_mon_to_interps = {};
                dec_interp_to_mons = {};
            end
            
            if detail_level > 1
                bounds = cell(size(prog.dec_polys));
                for i=1:size(prog.dec_polys, 1)
                    bounds{i} = repmat([-1, 1], prog.dec_polys(i).mon_basis.n, 1);
                end

                for i=1:size(prog.dec_polys, 1)
                    pol = prog.dec_polys(i);
                    dec_interp_basis = AffineFeketeBasis(pol.mon_basis.variables, ...
                                                        pol.mon_basis.d, bounds{i});
                    dec_interp_bases = [dec_interp_bases; dec_interp_basis];

                    [dec_interp_to_mon, dec_mon_to_interp] = dec_interp_basis.inter_to_mon(prog.dec_polys(i).mon_basis);
                    dec_mon_to_interps = [dec_mon_to_interps; dec_mon_to_interp];
                    dec_interp_to_mons = [dec_interp_to_mons; dec_interp_to_mon];

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
                As_int = cell(size(As));
                for i=1:size(As, 1)
                    for j=1:size(As, 2)
                        As_int{i, j} = -const_mon_to_interps{i} * As{i, j} * dec_interp_to_mons{j};
                    end
                end
            end
            
            fprintf("\nProgram Details:\n================\n\n")
            
            fprintf("Total number of indeterminates: %d\n-------------------------------\n", size(prog.indeterminates, 1))
            display(prog.indeterminates)
            
            fprintf("Total number of decision polynomials: %d\n-----------------------------------\n\n", size(prog.dec_bases, 1))
            
            for i=1:size(prog.dec_bases, 1)
                b = prog.dec_bases(i);
                fprintf("* Polynomial %d:\n\n", i)
                fprintf(" Variables:\n")
                display(b.variables)
                
                fprintf(" Degree: %d\n", b.d)
                
                fprintf(" Dimension: %d\n\n", b.num_monomials)
                
            end
            
            fprintf("Total number of SOS constraints: %d\n--------------------------------\n\n", size(prog.constraint_polys, 1))
            
            for i=1:size(prog.constraint_polys, 1)
                fprintf("* Constraint %d\n\n", i)
                fprintf(" Variables / bounds :\n\n")
                
                bounds = prog.K_wts{i};
                
                vars = prog.constraint_variables{i};
                
                names = name(vars);
                
                for j=1:size(vars, 1)
                    fprintf(" %s [%0.2f, %0.2f]\n", names{j}, bounds(j, 1), bounds(j, 2))
                end
                fprintf("\n SOS Cone Degree: %d\n\n", prog.constraint_degrees(i))
                
                if detail_level > 0
                    for j=1:size(As, 2)
                        fprintf(" size(A) for Decision Poly %d: %dx%d\n", j, size(As{i, j}))
                        fprintf(" Cond(A) for Decision Poly %d: %d\n", j, cond(As{i, j}))
                        if detail_level > 1
                            fprintf(" Cond(A_interp) for Decision Poly %d: %d\n\n", j, cond(As_int{i, j}))
                        end
                    end
                end
                
                fprintf(" Constraint A size: %dx%d\n\n", size(cell2mat(As_int(i, :))))
                
                fprintf(" Cond(A): %5d\n\n", cond(cell2mat(As_int(i, :))))
                
                
            end
            
            if detail_level > 1
                fprintf("\nsize(A): %dx%d\n", size(cell2mat(As_int)));
            end
            
            
            
        end
        
    end


end
