classdef FeketeBasis < InterpolantBasis & RectangularInterpolantBasis & handle & SerializableBasis
properties
    n double;
    d double;
    L double;
    U double;
    pts double;
    P0 double;
    P0_full;
    polynomials msspoly;
    variables msspoly;
    w double;
    bounds double;
end

properties (Constant)
    name = 'FeketeCube';
end

methods
    function int_params = FeketeBasis(n,d)
    % This method generates parameters for the interpolant basis representation
    % of sum-of-squares polynomials with an arbitrary number of variables.
    % --------------------------------------------------------------------------
    % USAGE of "FeketeCube"
    % intParams = FeketeCube(n,d)
    % --------------------------------------------------------------------------
    % INPUT
    % n:	number of arguments to the polynomials
    % d:	degree of polynomials to be squared
    %
    % OUTPUT
    % intParams:                interpolation parameters
    % - intParams.n:            number of arguments to the polynomials
    % - intParams.d:            degree of polynomials to be squared
    % - intParams.L:            dimension of the space of (intParams.n)-variate
    %                           degree-d/2 polynomials
    % - intParams.U:            dimension of the space of (intParams.n)-variate
    %                           degree-(d) polynomials
    % - intParams.pts:          approximate Fekete points for degree-(2)
    %                           polynomial interpolation. (intParams.U x 1) array.
    % - intParams.w:            (scaled) weights for Clenshaw-Curtis quadrature
    % - intParams.P0:           evaluations of bivariate product Chebyshev 
    %                           polynomials of the first kind up to degree d at
    %                           the points intParams.pts
    %                           (intParams.U x intParams.L) array.
    % - intParams.P:            evaluations of a basis for the space of
    %                           (intParams.n)-variate degree-d polynomials 
    %                           at the points intParams.pts. 
    %                           (intParams.U x intParams.L) array with 
    %                           orthonormal columns.  
    % --------------------------------------------------------------------------
    % EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
    % chebpolyval from Chebfun. Chebfun is an open-source package for  
    % numerical computing with functions: http://www.chebfun.org/.
    % Our code has been developed and tested with Chebfun version 5.5.0.
    % The latest version of the Chebfun package can be downloaded from
    % http://www.chebfun.org/download/.
    %
    % partitions. The partitions function computes all partitions of an integer.
    % We use the implementation of John D'Errico from a MATLAB Central File 
    % Exchange post which is available at
    % https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer.
    % -------------------------------------------------------------------------
        
        if check_for_basis(FeketeBasis.name, n, d)
            intParams = read_basis(FeketeBasis.name, n, d);
            int_params.populate_from_struct(intParams);
            return
        end

        intParams = FeketeBasis.generate_fekete_cube(n, d);
        'Writing to disk'
        write_basis(intParams);
        'Finished writing to disk'
        
        int_params.populate_from_struct(intParams)
        
    end
    
    
    function populate_from_struct(cube, dat)
        cube.n = dat.n;
        cube.L = dat.L;
        cube.U = dat.U;
        cube.pts = dat.pts;
        cube.P0_full = dat.P0_full;
        cube.P0 = dat.P0;
        cube.polynomials = dat.polynomials;
        cube.w = dat.w;

        if isfield(dat, 'mon_basis')
            cube.variables = dat.mon_basis.variables;
            dat.bounds = repmat([-1, 1], dat.mon_basis.n, 1);
            dat.d = dat.d * 2;
            cube.d = dat.d;
            cube.bounds = dat.bounds;
            dat.variables = dat.mon_basis.variables;
            dat = rmfield(dat, {'mon_basis', 'nrPoints', 'nrPoints1D', 'P', 'P_full'});
            overwrite_basis(dat);
        else
            cube.variables = dat.variables;
            cube.bounds = dat.bounds;
            cube.d = dat.d;

        end
        
    end

    function evaluations = eval_weight_polys(cube, pts)
        if size(pts, 2) ~= size(cube.vars)
            error('cube is %d dimensional, but pts is %d dimensional', size(pts, 2), size(cube.vars));
        end
        lb = -1 * ones(size(pts, 2), 1);
        ub = ones(size(pts, 2), 1);
        evaluations = bsxfun(@minus,pts,lb').*bsxfun(@minus,ub',pts);

    end

    function [P0_to_mon, mon_to_P0] = inter_to_mon(cube, mon_basis)
        if ~isequal(cube.variables, mon_basis.variables)
            error('AffineFeketeCube.variables does not match the provided monomial basis');
        end
        [P0_to_mon, mon_to_P0] = monomial_to_interpolant(cube.P0_full, cube.polynomials, mon_basis);
    
        % cube.P0_to_mon = P0_to_mon;
        % cube.P_to_mon = P_to_mon;
        % cube.mon_to_P0 = mon_to_P0;
        % cube.mon_to_P = mon_to_P;
    end


end

methods (Static)
    function intParams = generate_fekete_cube( n, d)
        if mod(d, 2) == 1
            error('Degree passed to FeketeCube() should be even');
        end
        name = 'FeketeCube';
        d = d/2;
        sprintf('Creating new %s with n=%d, d=%d', name, n, d)

        intParams.name = FeketeBasis.name;
        intParams.n = n;
        intParams.bounds = repmat([-1, 1], n, 1);
        intParams.d = 2*d;
        intParams.L = nchoosek(n+d,n);
        intParams.U = nchoosek(n+2*d,n);
        
        'creating monomial basis'

        mon_basis = MonomialBasis(intParams.n, intParams.d*2);
        intParams.variables = mon_basis.variables;

        'monomial basis generated'
        
        nrPoints1D = 2*d+1;

        'creating original point grid'
        nrPoints = prod(nrPoints1D:nrPoints1D+n-1);
        pts = zeros(nrPoints,n);
        for j = 1:n
            sprintf('starting x_%d', j)
            temp = 1;
            for i = 1:j-1; temp = kron(temp,ones(nrPoints1D+i-1,1)); end;
            temp = kron(temp,chebpts(nrPoints1D+j-1));
            for i = j+1:n; temp = kron(temp,ones(nrPoints1D+i-1,1)); end;
            pts(:,j) = temp;
        end

        'point grid generated'
        
        P = ones(nrPoints,intParams.U);
        m = ones(intParams.U,1);
        
        'evaluating product chebyshev polynomials'

        %prod_poly_vecs = 

        col     = 0;
        lrEye   = fliplr(eye(2*d+1)); %descending T_n


        mult_mats = cell(n, 2*d + 1);
        for var_ind = 1:n
            sprintf('generating variable %d', var_ind)
            for deg_ind=1:2*d + 1
                sprintf('generating degree %d', deg_ind)
                pol = chebyshev_in_variable(lrEye(:,deg_ind), ...
                                  var_ind, ...
                                  mon_basis);
                mult_mats{var_ind, deg_ind} = vector_poly_multiply_hack(pol, mon_basis, mon_basis);

            end
        end

        one_vec = msspoly_to_vector(msspoly(1), mon_basis);

        prod_polynomial_vecs = repmat(one_vec, 1, intParams.U);

        for t = 0:2*d      % polynomials with total degree up to 2*d
            sprintf('starting d = %d', t)
            allDegs = partitions(t, ones(1,n));
            [nrDegs,~] = size(allDegs);
            for i = 1:nrDegs
                sprintf('starting comb %d/%d', i, nrDegs)
                col = col+1;
                for j = 1:n
                    dj = allDegs(i,j);
                    prod_polynomial_vecs(:, col) = mult_mats{j, dj+1}  * prod_polynomial_vecs(:, col);
                    P(:,col) = P(:,col).*chebpolyval(lrEye(:,dj+1),pts(:,j));


                    if dj == 1; m(col) = 0;
                    else; m(col) = m(col)*(((-1)^dj+1)/(1-dj^2)); end;
                end
            end
        end
        'generating prod_polynomials'

        prod_polynomials = msspoly(ones(intParams.U, 1));
        for col_ind=1:intParams.U
            prod_polynomials(col_ind) = prod_polynomial_vecs(:, col_ind)' * mon_basis.monomials;
        end

        'Solving for weights'
        
        w = P'\m;
        ind = abs(w)>0;
        
        % extracts the positive entries of w
        w = w(ind);
        % extracts the subset of points indexed with the support of w
        pts = pts(ind,:);
        % extracts the subset of polynomials up to total degree d
        P0_large = P(ind, :);
        P = P(ind,1:intParams.L);

        'testing polynomials'
        %P' should equal P_test
        P_test = dmsubs(prod_polynomials, mon_basis.variables, pts');
        normm = norm(P_test - P0_large', 'fro');
        if norm(P_test - P0_large', 'fro') > 1E-8
            warning('Interpolant basis monomials may not be accurate\n ||P0 - P_test|| = %d, should be very small', normm);
        end


        intParams.w = w;
        intParams.pts = pts;
        intParams.P0 = P;

        'Doing QR nonsense'
        
        %[P_large, ~] = qr(P0_large);
        %[intParams.P,~] = qr(P,0);
        
        %intParams.P_full = P_large;
        intParams.P0_full = P0_large;
        

        
        intParams.polynomials = prod_polynomials;
        
        'Finished generating fekete cube'
    end

end

end