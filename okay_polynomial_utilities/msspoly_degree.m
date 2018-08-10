function degree = msspoly_degree(poly)
    [~, pows, ~] = decomp(poly);
    degree = max(sum(full(pows)'));
end