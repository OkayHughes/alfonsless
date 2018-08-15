base_path = '/home/owhughes/matlab';
addpath(fullfile(base_path, 'chebfun-master'));
addpath(fullfile(base_path, 'Padua2DM'));
addpath(fullfile(base_path, 'partitions'));
addpath(fullfile(base_path, 'spotless'));
addpath(fullfile(base_path, 'alfonso-master'));

addpath ./okay_poly_utils
addpath ./spotless_utils
addpath ./alfonso_utils
addpath ./interpolant_bases
addpath .
cwd = pwd;
cd(fullfile(base_path, 'spotless'));
spot_install_lite
cd(cwd);
