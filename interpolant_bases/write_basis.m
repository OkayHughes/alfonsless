function write_basis(basis)
    name = basis.name;
    n = basis.n;
    d = basis.d;
    ex = exist(name);
    if ex == 0
        mkdir("basis");
    elseif ex ~= 7
        error("%s is a non-directory", name);
    end
    ex = exist(fullfile(name, 'manifest.mat'));
    fname = sprintf("%d,%d.mat", n, d);
    write = false;
    if ex == 0
        manifest = [n, d];
        fnames = fname;
        save(fullfile(name, 'manifest.mat'), 'manifest', 'fnames');
        write = true;
    elseif ex ~= 2
        error("%s is the wrong format", fullfile(name, 'manifest.mat'));
    else
        rd = matfile(fullfile(name, 'manifest.mat'), 'writable', true);
        if ~ismember([n, d], rd.manifest, 'rows')
            rd.manifest = [rd.manifest;
                           [n, d]];
            rd.fnames = [rd.fnames;
                         fname];
            write = true;
        end
    end

    if write
        save(fname,'basis');
    end
end