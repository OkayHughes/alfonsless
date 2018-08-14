function basis = read_basis(name, num_vars, degree)
    if ~check_for_basis(name, num_vars, degree)
        error("Basis with name %s has not been generated for n=%d, d=%d", name, num_vars, degree);
    end

    mname = fullfile(name, 'manifest.mat');
    rd = matfile(mname);
    manifest = rd.manifest;
    [~, idx] = ismember([n, d], manifest, 'rows');
    fname = rd.fnames(idx, 1);
    obj = matfile(fullfile(name, fname));
    basis = obj.basis;
end