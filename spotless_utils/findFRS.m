function out = findFRS(prob)
%% extract parameter
    disp('Extracting parameters')
    t = prob.t ;
    z = prob.z ;
    k = prob.k ;
    f = prob.f ;
    hZ = prob.hZ ;
    hK = prob.hK ;
    hZ0 = prob.hZ0 ;
    degree = prob.degree ;
    domain_size = prob.domain_size ;

    if isfield(prob,'costs')
        costs = prob.costs ;
    elseif isfield(prob,'int_ZK')
        int_ZK = prob.int_ZK ;
    end

    if isfield(prob, 'g')
        g = prob.g ;
    end

    if isfield(prob, 'z2track')
        z2track = prob.z2track ;
    else
        z2track = z ;
    end

    if isfield(prob, 'hZ2track')
        hZ2track = prob.hZ2track ;
    else
        hZ2track = hZ ;
    end

    if isfield(prob,'hDelta')
        hDelta = prob.hDelta ;
    end

    if isfield(prob,'run_solver')
        run_solver = prob.run_solver ;
    else
        run_solver = 0 ;
    end
    
    if isfield(prob,'enforce_boundary')
        enforce_boundary = prob.enforce_boundary ;
    else
        enforce_boundary = 0 ;
    end
    
    if isfield(prob,'reduce_vmon')
        reduce_vmon = prob.reduce_vmon ;
    else
        reduce_vmon = 0 ;
    end

%% define variables
    disp('Defining problem variables')

    % initialize program and indeterminate vars
    prog = spotsosprog;
    prog = prog.withIndeterminate(t); % Time
    prog = prog.withIndeterminate(z); % States
    prog = prog.withIndeterminate(k); % Params

    % create v polynomial
    if reduce_vmon
        vmontz = monomials([t;z], 0:degree) ;
        vmonzk = monomials([z;k], 0:degree) ;
        vmontk = monomials([t;k], 0:degree) ;
        vmon = [vmontz; vmonzk; vmontk] ;
        vmonspec = {[t;z], [z;k], [t;k]} ;
    else
        vmon = monomials([t;z;k], 0:degree) ;
    end
    [prog, v, ~] = prog.newFreePoly(vmon) ;

    % create w polynomial
    wmon = monomials([z2track;k], 0:degree) ;
    [prog, w, wcoeff] = prog.newFreePoly(wmon) ;

    % if disturbance is included, create q
    if exist('g','var')
        q = msspoly(zeros([size(g,2),1])) ;
        qmon = monomials([t;z;k], 0:degree) ;
        for qidx = 1:length(q)
            [prog, q(qidx), ~] = prog.newFreePoly(qmon) ;
        end
    end

    % create variables for the equality constraints of the dual formulation
    v0 = subs(v,t,0) ;
    dvdt = diff(v,t) ;
    dvdz = diff(v,z) ;
    Lfv = dvdt + dvdz*f ;
    if exist('g','var')
        Lgv = dvdz*g ;
    end

    % plug in delta conditions to v, which will capture the system trajs,
    % and create zidx which indexes all the non-delta initial states
    zidx = 1:length(z) ;
    if exist('hDelta','var')
        [nDelta,~] = size(hDelta) ;
        for idx = 1:nDelta
            v0 = msubs(v0,z(hDelta(idx,1)), hDelta(idx,2)) ;
        end
        zidx(hDelta(:,1)) = [] ;
    end

%% define constraints
    disp('Defining constraints')
    tic

    % define the time range T = [0,1]
    hT = t*(1-t);

    if exist('g','var')
        % -Lfv - q > 0 on T x Z x K
        prog = sosOnK(prog, -Lfv - ones(size(q))'*q, [t;z;k], [hT; hZ; hK], degree) ;
    else
        % -Lfv > 0 on T x Z x K
        prog = sosOnK(prog, -Lfv, [t;z;k], [hT; hZ; hK], degree) ;
    end

    % v(t,.) + w > 1 on T x Z x K
    if reduce_vmon
        prog = sosOnK(prog, v + w - 1, [t;z;k], [hT; hZ; hK], degree, vmonspec) ;
    else
        prog = sosOnK(prog, v + w - 1, [t;z;k], [hT; hZ; hK], degree) ;
    end

    % w > 0 on Z x K
    prog = sosOnK(prog, w, [z2track;k], [hZ2track; hK], degree) ;

    % -v(0,.) > 0 on Z0 x K
    prog = sosOnK(prog, -v0, [z(zidx);k], [hZ0; hK], degree) ;

    % v(t,.) > 0 on T x dX x K if boundary is active
    if enforce_boundary
        for dzidx = 1:length(hZ)
            prog = sosOnK(prog, v, [t;z;k], [hT; hZ ; -hZ(dzidx); hK], degree) ;
        end
    end
    
    % if disturbance is included, we need the following constraints:
    % q - Lgv > 0 on T x Z x K
    % q + Lgv > 0 on T x Z x K
    % q > 0 on T x Z x K
    if exist('g','var')
        for qidx = 1:length(q)
            prog = sosOnK(prog, q(qidx) - Lgv(qidx), [t;z;k], [hT;hZ;hK], degree) ;
            prog = sosOnK(prog, q(qidx) + Lgv(qidx), [t;z;k], [hT;hZ;hK], degree) ;
            prog = sosOnK(prog, q(qidx), [t;z;k], [hT;hZ;hK], degree) ;
        end
    end    

    toc

%% generate cost function
    disp('Generating cost function')

    tic
    if exist('int_ZK','var')
        obj = int_ZK(wmon)' * wcoeff ;
    elseif exist('costs','var')
        % The variable 'costs' has the following format:
        % {vars1, cost_func1(vars1) ;
        %  vars2, cost_func2(vars2) ;
        %          ...
        %  varsn, cost_funcn(varsn)}

        % initialize total cost function as 0
        obj = 0 ;

        [n_cfuncs,~] = size(costs) ;

        if n_cfuncs > 1
            % if there is more than one cost function to mix
            % for each monomial...
            for widx = 1:length(wmon)
                % initialize the product of the cost functions as for the
                % current monomial
                objmon = 1 ;

                % for each cost function and its variables...
                for cidx = 1:n_cfuncs
                    % extract the particular variables and cost function
                    cvars = costs{cidx,1} ;
                    cfunc = costs{cidx,2} ;

                    % get the complement of the variables corresponding to the
                    % current cost function; this is needed because the second
                    % argument to 'decomp' is what variables must be *excluded*
                    % from the decomposition
                    ztemp = decomp([z;k], cvars) ;

                    % check if the current vars are in the current monomial by
                    % excluding their complement
                    [cmon,pmon,~] = decomp(wmon(widx), ztemp) ;

                    % if the variables are in the current monomial...
                    if ~isempty(cmon)
                        % recompose the current variables with their powers
                        ccomp = recomp(cmon,pmon,1) ;

                        % evaluate the cost function with the relevant
                        % chunk of the monomial
                        costmon = cfunc(ccomp) ;
                    else
                        % ... or just leave the cost contribution as 1
                        costmon = 1 ;
                    end

                    % multiply the current monomial's cost function by
                    % the current cost
                    objmon = objmon * costmon ;
                end

                % add the current monomial's cost to the total obj function
                obj = obj + wcoeff(widx)*objmon ;
            end
        else
            cfunc = costs{1,2} ;
            obj = cfunc(wmon)' * wcoeff ;
        end
    end
    toc

    out.prog = prog ;
    out.wmon = wmon ;
    out.wcoeff = wcoeff ;
    out.obj_func = obj ;

    if run_solver
    %% set options
        options = spot_sdp_default_options();
        options.verbose = 1;
        options.domain_size = domain_size;
        options.solveroptions = [];

    %% solve robustness problem
        disp('Solving the robustness problem')
        tic
        if strcmp(prob.solver,'mosek')
            sol = prog.minimize(obj, @spot_mosek, options);
        elseif strcmp(prob.solver,'scs')
            sol = prog.minimize(obj, @spot_scs, options);
        elseif strcmp(prob.solver,'sdpt3')
            options.sdpt3_pinfeas_tol = 1e-5 ;
            options.sdpt3_dinfeas_tol = 1e-5 ;
            sol = prog.minimize(obj, @spot_sdpt3, options);
        elseif strcmp(prob.solver,'sedumi')
            sol = prog.minimize(obj, @spot_sedumi, options);
        elseif strcmp(prob.solver,'frlib')
            options.approx = 'dd' ;
            sol = prog.minimize(obj, @spot_frlib, options);
        end
        end_time = toc ;

    %% extract problem info
        disp('Extracting results')
        out.v = sol.eval(v) ;
        out.w = sol.eval(w) ;
        out.final_cost = sol.eval(obj) ;
        out.duration = end_time ;
        if exist('g','var')
            out.qs = sol.eval(q);
        end

        disp([num2str(end_time/3600), ' hrs elapsed solving problem'])
    end
end
