function basis = monomial_basis(num_variables, max_degree, variables)
    %monomial_basis(num_variables, max_degree, var_array)
    %Takes:
    %   num_variables: integer
    %       number of variables over which the polynomials can
    %       be expressed
    %   max_degree: integer
    %       the maximum degree of polynomials in the
    %       space we're generating a basis for
    %   variables=var_name:
    %       string to be used as the variable name
    %   ||variables=var_array: nx1 msspoly vector
    %       msspoly array to be assigned as basis.variables
    %Returns:
    %   basis: struct
    %       basis.n: integer
    %           number of variables in the polynomial
    %       basis.d:
    %           highest degree of the monomial
    %       basis.num_monomials
    %           nchoosek(n + d, d)
    %       basis.variables: mss
    %           num_variables x 1 vector of free msspoly
    %           variables, the output of msspoly(var_name, num_variables)
    %       basis.monomials: msspoly variable vector
    %           nchoosek(max_degree + num_variables, max_degree)
    %           x 1 vector of msspoly monomials, the output of 
    %           monomials(basis.variables, 0:max_degree)
    %       basis.power_matrix: real matrix
    %           nchoosek(max_degree + num_variables, max_degree) x n matrix
    %           A, the output of monomials_to_matrix(basis.monomials,
    %           basis.variables).
    basis.n = num_variables;
    basis.d = max_degree;
    if isa(variables, 'char')
        basis.variables = msspoly(variables, num_variables);
    elseif isa(variables, 'msspoly')
        if size(variables, 1) ~= basis.n
            error('array passed as `variables` has size %d, but n=%d', size(variables, 1), basis.n)
        end
        basis.variables = variables;
    end
    basis.monomials = monomials(basis.variables, 0:max_degree);
    basis.num_monomials = size(basis.monomials, 1);
    basis.power_matrix = monomials_to_matrix(basis.monomials, basis.variables);
end