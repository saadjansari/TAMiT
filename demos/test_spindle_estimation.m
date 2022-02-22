% DEMO: Spindle Line Detection
fprintf('#########    DEMO : Spindle Line Detection    ##########\n\n')

% Add relative paths
addpath('../classes');


%% Load intensity image

% Path to .tif file of spindle intensity image
image_path = './image_spindle.tif';

% Read the tif stack
img3D = read_tiff_stack(image_path);

% To view this image, uncomment, and run the following lines
% figure; imagesc( max(img3D,[],3)); 
% title('Intensity image (max Z projection')

%% Define params

% Region Finding: Expected Minor Axis Length
params.expectedMAL = 5;

% Region Finding: Minimum area (in pixels)
params.minRegionArea = 6;

% Thickness of spindle axis used to calculate spindle intensity
params.linewidth = 3;

% A setting (either 0 or 1) that forces one of the ends of the spindle to
% be at the brightest pixel. 
params.brightestPixelAsSPB = 0;

% Other parameters
params.visuals = 1; % Visuals ON. Necessary to produce demo figures.
params.visuals_path = './demo_figs'; % Path of folder where to save the demo figures.
params.verbose = 1; % Command line outputs maximized.
% -------------------

%% Run spindle finder

% Delete visuals path before starting. This ensures a clean directory.
if exist(params.visuals_path) == 7
    rmdir(params.visuals_path,'s')
end

% Run spindle detection function used in TAMiT. 
MitoticCell.findTheSpindle(img3D, params);
