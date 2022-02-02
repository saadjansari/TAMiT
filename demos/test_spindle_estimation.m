% TEST
% Spindle Line Detection


fprintf('#########    Test: Spindle Line Detection    ##########\n\n')

% Load Mitotic Cell SPB image
img3D = read_tiff_stack('image_spindle.tif');

% Add relative paths
addpath('../classes');

% Define params
params.visuals = 1;
params.visuals_path = './demo_figs';
params.verbose = 1;
params.linewidth = 3;
params.brightestPixelAsSPB = 0;
params.spindleDeterminationSensitivity = 0.4;
params.spindleMinIntensity = 0.35;

% Delete visuals path
if exist(params.visuals_path) == 7
    rmdir(params.visuals_path,'s')
end

% Run detection
MitoticCell.findTheSpindle(img3D, params);
