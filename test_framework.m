function res = test_framework(deg, num_consts, defectiveness)
	addpath('okay_sample_probs');
	setup;

	consts = cell(num_consts, 3);
    
	val = nchoosek(deg + 5, 5);
	for i=1:4
		for j=1:num_consts
			if (defectiveness == 0)
				A = rand(val);
			elseif (defectiveness == 1)
				A = rand(val);
				A(:, end) = 0;
			else
				A = zeros(val);
				for i=1:(floor(val/2))
					A(i*2, :) = rand(val, 1);
				end
			end
			consts{j, i} = A;
		end
    end
    res = SOS_fix_deg_exten(deg, consts);
end


