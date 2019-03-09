function out = call_exp(degree, indices, name)
    addpath segway
    setup
    knockout = zeros(7,1);
    knockout(indices) = 1;
    diary(name)
    out = segway_FRS_solver_okay(degree, knockout);
    diary off
end