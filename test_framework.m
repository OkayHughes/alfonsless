function res = test_framework(num_polys, deg, num_consts)
	addpath('okay_sample_probs');
	setup;

	consts = cell(num_consts, num_polys);
    
	val = nchoosek(deg + 5, 5);
	for i=1:num_polys
		for j=1:num_consts
			A = rand(val);
			consts{j, i} = A;
		end
    end
    res = SOS_fix_deg_exten(num_polys, deg, consts);
end


