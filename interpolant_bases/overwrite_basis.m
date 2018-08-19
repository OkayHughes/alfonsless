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
    if ex ~= 2
        error('%s is the wrong format', fullfile(name, 'manifest.mat'));
    else
        rd = matfile(fullfile(prefix, name, 'manifest.mat'), 'writable', true);
        if ismember([n, d], rd.manifest, 'rows')
             write = true;
        else
          warning('File to be overwritten did not exist')
        end
    end

    if write
      try
        class(basis)
        save(fullfile(prefix, name, fname),'basis', '-v7.3');
      catch ME
        warning('%s with n=%d, d=%d', name, n, d);
        return 
      end

    end
end
