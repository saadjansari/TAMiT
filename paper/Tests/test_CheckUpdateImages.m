% Test Script for CheckUpdateImages functionality
clc; close all;

opts.CFG = 'RELEASE';
opts.LOC = 'Local';

warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath( genpath(pwd) );
warning('on', 'MATLAB:rmpath:DirNotFound');
addpath( pwd);

% check if params file exists in the current folder
if exist( fullfile( pwd, 'initParams.m' ) ) ~= 2
    error( ['copy initParams.m to the current folder location : ', pwd] );
end

% run the settings file : creates a params.mat file in the save directory folder
paramsPath = feval('initParams', opts);

% Define cleanup tasks
c2 = onCleanup( @() delete( gcp('nocreate') ) );
c3 = onCleanup( @() disp('Closing files and cleaning up') );

try; paramsPath = paramsPath{1}; end

addpath( genpath( [pwd, filesep, 'functions'] ) );
addpath( [pwd, filesep, 'classes'] );

% Load the params
load( paramsPath); 

% Create fake movie
img = zeros(150,150,7,25,2);
img(:,:,:,8,2) = randn(150,150,7);
img(:,:,:,15,2) = randn(150,150,7);
img(:,:,:,22,2) = randn(150,150,7);
imageData = ImageData( img, 'Lifetime', [5 25]); 

% Run the specific type of cell
switch params.cellInfo.type
    case 'Monopolar'

        % Initialize Monopolar Cell
        myCell = MonopolarCell( imageData, ...
            params.cellInfo.channelFeatures, ...
            params.cellInfo.channelsToFit, ...
            params, ... 
            'Species', params.cellInfo.species, ...
            'Strain', params.cellInfo.strain);

        [myCell, imgNew] = CheckUpdateImages( myCell, 2);

        % play movie of old 
        imageData.Play()
        pause(1)

        % play movie of new 
        imageDataNew = ImageData( imgNew, 'Lifetime', [5, 25]); 
        imageDataNew.Play()

end


