function res = check_for_basis(name, n, d)
    prefix = fileparts(which(mfilename));
    ex = exist(fullfile(prefix, name));
    if ex ~= 7
        res = false;
        return;
    end

    ex = exist(fullfile(prefix, name, 'manifest.mat'));

    if ex ~= 2
        res = false;
        return;
    else
        rd = matfile(fullfile(prefix, name, 'manifest.mat'));
        if ismember([n, d], rd.manifest, 'rows')
            res = true;
            return;
        else
            res = false;
            return;
        end
    end
end