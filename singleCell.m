function features = singleCell( paramsPath)
    % Input: paramsPath is a path to a .mat file

    clc; close all;
    clearvars -except paramsPath 

    % Make sure we have the correct paths
%     warning('off', 'MATLAB:rmpath:DirNotFound');
%     rmpath( genpath(pwd) );
%     warning('on', 'MATLAB:rmpath:DirNotFound');
    addpath( pwd);
    addpath( genpath( [pwd, filesep, 'functions'] ) );
    addpath( [pwd, filesep, 'classes'] );

    % Load the params
    load( paramsPath); 

    % Start diary
    c1 = onCleanup( @() eval('diary off') );
    diary( fullfile( params.saveDirectory, 'singleCell.log') );

    disp('--------------------------------------------------------------------------------------')
    disp( ['Configuration: ' params.CFG])
    disp(' ')
    disp('Paths:')
    disp(' ')
    disp( ['    Run Path: ' params.runPath])
    disp( ['    Save Path: ' params.saveDirectory])
    disp( ['    User Settings: ' params.paramsPath, '.mat' ])
    disp(' ')
    disp('--------------------------------------------------------------------------------------')

    % Locate the single cell movie
    if isempty( params.cellInfo.moviePath) 
        [cfile, cpath] = uigetfile({'*.mat';'*.tiff'});
        if cfile == 0
            error( 'singleCell: no file was selected for import'); end
        params.moviePath = [cpath, cfile];
    end

    % Initialize image data from single cell movie
    if ~isempty(params.cellInfo.lifetime)
        imageData = ImageData.InitializeFromCell( params.cellInfo.moviePath, ...
            'Lifetime', params.cellInfo.lifetime );
    else
        imageData = ImageData.InitializeFromCell( params.cellInfo.moviePath);
    end
    sizeVoxels = imageData.GetSizeVoxels;
    timeStep = imageData.GetTimeStep;
    save( paramsPath, 'sizeVoxels', 'timeStep', '-append')

    % Run the specific type of cell
    switch params.cellInfo.type
        case 'Mitosis'
            fname = 'MitoticCell';
        case 'Monopolar'
            fname = 'MonopolarCell';
        case 'Interphase'
            fname = 'InterphaseCell';
    end
    
    myCell = feval( fname, imageData, ...
        params.cellInfo.channelFeatures, ...
        params.cellInfo.channelsToFit, ...
        params, ... 
        'Species', params.cellInfo.species, ...
        'Strain', params.cellInfo.strain);
   
            
    % Find features
    myCell = myCell.FindFeatures();


end


