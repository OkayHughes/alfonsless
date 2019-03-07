for I in 2 4
do
	for J in 1 5 10
	do
		for K in 4 8 10 
		do
			./test_grid.py $I $K $J "3-3"&
		done
	done
done
