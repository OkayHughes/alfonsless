function res = test_segway(deg)
	addpath('segway');
	setup;
	
	res = segway_FRS_solver_okay(deg);
end


