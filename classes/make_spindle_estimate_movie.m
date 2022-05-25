% TEST: SPINDLE DETECTION

%% PARAMS

% Path to segmented movie
mov_path = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/LukeTest/1149_100R_100G_25_deg_001_A.mat';

% Specify channel to analyze
mt_channel = 1;

% Estimation params

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
% Set path to save figs
params.visuals_path = './test_figs'; % Path of folder where to save the demo figures.
params.verbose = 0; % Command line outputs maximized.
% -------------------

%% ESTIMATION 

fprintf('#########    TEST : Spindle Line Detection    ##########\n\n')

% Add relative paths
addpath(genpath('../../classes/.'));

% Load movie
imgData = ImageData.InitializeFromCell( mov_path);
imgData.PrintInfo()

% Play Movie (optional)
% imgData.Play()

% Extract images for MT channel
imXYZTC = imgData.GetImage();
if length( size( imXYZTC) ) == 5
    imXYZT = imXYZTC(:,:,:,:,mt_channel);
elseif length( size( imXYZTC) ) == 4
    imXYZT = imXYZTC(:,:,:,:);
else
    error('dimensionality is unexpected')
end
clear imXYZTC

% Delete visuals path before starting. This ensures a clean directory.
if exist(params.visuals_path) == 7
    rmdir(params.visuals_path,'s')
end

%% Estimate Spindle

nFrames = size(imXYZT,4);

% Set up movie making
vidfile = VideoWriter( 'spindle_estimate.mp4','MPEG-4');
vidfile.FrameRate = 5;
open(vidfile);

% Loop over time
for jFrame = 1: nFrames
    
    disp(['Frame = ', num2str(jFrame)])
    if exist(params.visuals_path) == 7
        rmdir(params.visuals_path,'s')
    end
    
    % 3D image
    img3D = imXYZT(:,:,:,jFrame);
    
    % Estimate Spindle
    % Run spindle detection function used in TAMiT. 
    [spindle,h] = MitoticCell.findTheSpindle(img3D, params);
    sgtitle(['Frame = ',num2str(jFrame)]) 
    writeVideo(vidfile, getframe(h));
    close(h)
    
end
close(vidfile)
