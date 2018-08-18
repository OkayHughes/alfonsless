classdef (Abstract) InterpolantBasis

properties(Abstract)
    name char;
    n double;
    d double;
    L double;
    U double;
    pts double;
    P0 double;
    polynomials msspoly;
    variables msspoly;
end

methods(Abstract)
    eval_weight_polys(prog, pts);
    inter_to_mon(prog, mon_basis);
end

end