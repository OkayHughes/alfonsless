cwd = pwd;
try
    base_path = 'Path to directory containing dependencies';
    addpath(fullfile(base_path, 'chebfun-master'));
    addpath(fullfile(base_path, 'Padua2DM'));
    addpath(fullfile(base_path, 'partitions'));
    addpath(fullfile(base_path, 'spotless'));
    addpath(fullfile(base_path, 'alfonso-master')); 

    addpath ./poly_utils
    addpath ./spotless_utils
    addpath ./alfonso_utils
    addpath ./interpolant_bases
    addpath ./alfonsless
    addpath .
    cd(fullfile(base_path, 'spotless'));
    spot_install_lite
    cd(cwd);
catch e
   fprintf(1,'%s\n',e.identifier);
   fprintf(1,'%s\n',e.message);
   cd(cwd);
end