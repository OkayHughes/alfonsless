degree = 12 ;
intParams = FeketeCube(1,degree) ;
basis = FeketeCube.mon_basis
% point constraint in the problem description
p = 0; 
c = 1.5 ;

%% SPOTLESS PROBLEM
% polynomial variable
t = basis.variables;

for i=1:100
    vec = rand(size(basis.monomials));
    P0 = intParams.mon_to_P0 * vec;
    mon = vec' * basis.monomials;
    vals = msubs(mon, t, intParams.pts')
    norm(P0 - vals)
end