function difference = var_set_minus(minuend, subtrahend)
    isnt_in = @(variable) var_to_index(variable, subtrahend) == 0;
    var_inds = arrayfun(isnt_in, minuend);
    difference = minuend(var_inds);

end