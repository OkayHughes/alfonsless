function basis = read_basis(name, n, d)
    prefix = fileparts(which(mfilename));
    if ~check_for_basis(name, n, d)
        error('Basis with name %s has not been generated for n=%d, d=%d', name, n, d);
    end

    mname = fullfile(prefix, name, 'manifest.mat');
    rd = matfile(mname);
    manifest = rd.manifest;
    [~, idx] = ismember([n, d], manifest, 'rows');
    fnames = rd.fnames;
    obj = matfile(fullfile(prefix, name, char(fnames{idx, 1})));
    basis = obj.basis;
end