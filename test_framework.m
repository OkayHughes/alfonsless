function res = test_framework(deg, num_consts)
	addpath('okay_sample_probs');
	setup;

	consts = cell(num_consts, 3);
    
	val = nchoosek(deg + 5, 5);
	for i=1:4
		for j=1:num_consts
			A = rand(1, val);
			consts{j, i} = A;
		end
    end
    res = SOS_fix_deg_exten(deg, consts);
end


