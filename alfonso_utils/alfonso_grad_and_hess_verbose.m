function [in, g, H, L] = alfonso_grad_and_hess_verbose(x, params)
% This method computes the gradient and Hessian of the barrier function for
% the problem of polynomial envelopes.
% --------------------------------------------------------------------------
% USAGE of "alfonso_grad_and_hess"
% [in, g, H, L] = alfonso_grad_and_hess(x, params)
% --------------------------------------------------------------------------
% INPUT
% x:                    primal iterate
% params:               parameters associated with the method gH
% - params.numPolys:    number of decision polynomials
% - params.bnu:         complexity parameter of the augmented barrier (nu-bar)
% - params.n_arr:           number of arguments to the polynomials
% - params.d_arr:           largest degree of polynomials to be squared
% - params.U_arr:           dimension of the space of (params.n)-variate 
%                       degree-(2*params.d) polynomials

% - params.L_arr:           dimension of the space of (params.n)-variate 
%                       degree-(params.d) polynomials.
% - params.LWts_cell:        dimensions of the "weighted" polynomial spaces. 
%                       params.LWts(j) is the dimension of the space of
%                       (params.n)-variate degree-(params.d-1) polynomials
%                       for j = 1,...,n.
% - params.P_cell:           evaluations of the basis for the space of 
%                       (params.n)-variate degree-(params.d) polynomials
%                       at the interpolation points
% - params.PWts_cell:        evaluations of "weighted" basis polynomials at the
%                       interpolation points. params.PWts{j} is
%                       the evaluations of the basis   for the "weighted"
%                       space of (params.n)-variate degree-(params.d-1)
%                       polynomials at the interpolation points for 
%                       j = 1,...,n. the weight corresponding to 
%                       params.PWts{j} is sqrt(1-t_j^2).
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the 
%       interior of the cone.
% g:    gradient of the barrier function at x
% H:    Hessian of the barrier function at x
% L:    Cholesky factor of H
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------


    fprintf("\nENTERING GRADIENT AND HESSIAN CALCULATION\n");
    fprintf("=========================================\n\n");
    
    U_arr = params.U_arr;
    in = 1;
    g = zeros(sum(U_arr),1);
    H = zeros(sum(U_arr)); 
    L = zeros(sum(U_arr));
    
    numPolys = params.numPolys;

    % ORDER OF VARIABLES 
    % x = [x_1; x_2; ...] \in WSOS_(n,2*d)^* x WSOS_(n,2*d)^* x ...
    % x_1 corresponds to the 1st approximated polynomial, x_2 to the 2nd,...
    off = 0;
    for polyId = 1:numPolys
        fprintf("Decision Polynomial %d:\n", polyId);
        fprintf("-----------------------\n\n");
        n = params.n_arr(polyId);
        U = params.U_arr(polyId);
        P = params.P_cell{polyId};
        PWts = params.PWts_cell{polyId};
        xPoly = x(off+(1:U));
        % for the weight 1
        fprintf("H_1(x):\n");
        fprintf("--------------------\n\n");
        [inPoly, gPoly, HPoly] = gH_SOSWt(xPoly,P);

        if inPoly == 1
            for j = 1:n
                fprintf("H_%d(x):\n", j);
                fprintf("--------------------\n\n");
                % for the weight 1-t_j^2
                [inPolyWt, gPolyWt, HPolyWt] = gH_SOSWt(xPoly,PWts{j});
                inPoly  = inPoly & inPolyWt;
                if inPoly == 1
                    gPoly   = gPoly+gPolyWt;
                    HPoly   = HPoly+HPolyWt;
                    fprintf("Cond(H_{1..%d}(x)):  %d\n\n", j, cond(HPoly));
                else
                    gPoly   = NaN;
                    HPoly   = NaN;
                    break;
                end
            end
        end

        if inPoly == 1
            % checks positive semidefiniteness of HPoly one more time.
            % HPoly may fail Cholesky factorization even if inPoly == 1
            % due to numerical errors in summing up HPolyWt's above.
            [LPoly, err] = chol(HPoly, 'lower');
            inPoly  = inPoly & (err == 0);
        end

        if inPoly == 1
            g(off+(1:U)) = gPoly;
            H(off+(1:U),off+(1:U)) = HPoly;
            L(off+(1:U),off+(1:U)) = LPoly;
            off = off + U;
        else
            in = 0;
            g = NaN;
            H = NaN;
            L = NaN;
            warning("Cholesky factorization failed for a polynomial!");
            return;
        end
    end
    fprintf("Cond(H(x)) total: %5d\n", cond(H));
    fprintf("Cond(L) total: %5d\n\n", cond(L));
end

function [in, g, H] = gH_SOSWt(x, P)
% This method computes the gradient and Hessian of a barrier function term
% corresponding to a single weight and single approximated polynomial for
% the problem of polynomial envelopes.
% --------------------------------------------------------------------------
% USAGE of "gH_SOSWt"
% [in, g, H] = gH_SOSWt(x, P)
% --------------------------------------------------------------------------
% INPUT
% x:    subvector of the primal iterate corresponding to a single
%       approximated polynomial
% P:    evaluations of "weighted" basis polynomials at the
%       interpolation points
%
% OUTPUT
% in:   0 if P'*diag(x)*P is not positive definite. 1 if P'*diag(x)*P is
%       positive definite.
% g:    gradient of the barrier function term corresponding to a single
%       weight and single approximated polynomial at x
% H:    Hessian of the barrier function term corresponding to a single
%       weight and single approximated polynomial at x
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    Y = P'*diag(x)*P;
    fprintf("F(x): %5d\n", log(det(Y)));
    if ~issymmetric(Y)
        Y = (Y+Y')/2;
    end
    fprintf("Cond(Lambda(x)): %5d\n", cond(Y));
    [L, err] = chol(Y, 'lower');
    if err > 0
        in = 0;
        g = NaN;
        H = NaN;
    else
        in = 1;
        fprintf("Cond(L): %5d\n", cond(L));
        V = L\P';
        VtV = V'*V;

        g = -diag(VtV);
        H = VtV.^2;
        fprintf("Cond(H_i(x)): %5d\n", cond(H));
    end
end