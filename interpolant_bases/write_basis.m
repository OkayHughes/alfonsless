function write_basis(basis)
    prefix = fileparts(which(mfilename));
    name = basis.name;
    n = basis.n;
    d = basis.d;
    ex = exist(fullfile(prefix, name));
    if ex == 0
        mkdir(name);
    elseif ex ~= 7
        error('%s is a non-directory', name);
    end
    ex = exist(fullfile(prefix, name, 'manifest.mat'));
    fname = sprintf('%d,%d.mat', n, d);
    write = false;
    if ex == 0
        manifest = [n, d];
        fnames = {string(fname)};
        rd = matfile(fullfile(prefix, name, 'manifest.mat'), 'writable', true);
        write = true;
    elseif ex ~= 2
        error('%s is the wrong format', fullfile(name, 'manifest.mat'));
    else
        rd = matfile(fullfile(prefix, name, 'manifest.mat'), 'writable', true);
        if ~ismember([n, d], rd.manifest, 'rows')
             fnames = [rd.fnames; string(fname)];
             manifest = [rd.manifest; [n, d]];
             write = true;
        end
    end
    
    if write
      try
        save(fullfile(prefix, name, fname),'basis', '-v7.3');
      catch ME
        warning('%s with n=%d, d=%d', name, n, d);
        return 
      end
      rd.manifest = manifest;
      rd.fnames = fnames;

    end
end
