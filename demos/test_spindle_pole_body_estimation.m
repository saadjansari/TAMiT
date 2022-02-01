% TEST
% Spindle Poly Body Detection


fprintf('#########    Test: Spindle Poly Body Detection    ##########\n\n')

% Load Mitotic Cell SPB image
img3D = read_tiff_stack('image_spb.tif');

% Add relative paths
addpath('../classes');

% Define params
params.visuals = 1;
params.visuals_path = './demo_figs';
params.verbose = 1;
params.min_spot_area = 25;
params.min_conn_region_area = 7;

% Delete visuals path
if exist(params.visuals_path) == 7
    rmdir(params.visuals_path,'s')
end

% Run detection
MitoticCell.findSPB_sid4( img3D, params);