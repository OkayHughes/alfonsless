function gH_Params = gen_grad_params(intParams_arr)

    numPolys = size(intParams_arr, 2);

    gH_Params.numPolys = numPolys;
    gH_Params.n_arr = zeros(numPolys, 1);
    gH_Params.d_arr = zeros(numPolys, 1);
    gH_Params.U_arr = zeros(numPolys, 1);
    gH_Params.L_arr = zeros(numPolys, 1);
    gH_Params.LWts_cell = cell(numPolys, 1);
    gH_Params.P_cell = cell(numPolys, 1);
    gH_Params.PWts_cell = cell(numPolys, 1);
    gH_Params.bnu = 0;

    for i=1:numPolys
        intParams_i = intParams_arr(i);
        n   = intParams_i.n;
        d   = intParams_i.d;
        U   = intParams_i.U;
        L   = intParams_i.L;
        P   = intParams_i.P0;
        pts = intParams_i.pts;
        bounds = intParams_i.bounds;
        lb = bounds(:, 1);
        ub = bounds(:, 2);

        % dimension of the weight polynomial space (should be dimension of d)
        LWts = repmat(nchoosek(n+d-1,n),n,1);
        % parameter object for cone gradient/hessian function
        gH_Params.n_arr(i)      = n;
        gH_Params.d_arr(i)      = d;
        gH_Params.U_arr(i)      = U;
        gH_Params.L_arr(i)      = L;
        gH_Params.LWts_cell{i}  = LWts;
        nu                      = L+sum(LWts) ;
        gH_Params.bnu           = gH_Params.bnu + nu;
        gH_Params.P_cell{i}     = P;

        % create polynomial hT (g in the alfonso paper) to define space T = [-1,1]^2

        wtVals  = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);

        PWts = cell(n,1);
        for j = 1:n
            PWts{j}         = diag(sqrt(wtVals(:,j)))*P(:,1:LWts(j));
            [PWts{j}, ~]    = qr(PWts{j}, 0);
            % associated positive semidefinite cone constraints: 
            % PWts{j}'*diag(x_1)*PWts{j} >= 0,
            % PWts{j}'*diag(x_2)*PWts{j} >= 0,...
        end
        gH_Params.PWts_cell{i} = PWts;
    end
    gH_Params.bnu = gH_Params.bnu + 1;

end