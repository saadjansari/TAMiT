function features = singleCell( paramsPath)
    % Input: paramsPath is a path to a .mat file

    clc; close all;
    clearvars -except paramsPath 

    % Make sure we have the correct paths
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath( genpath(pwd) );
    warning('on', 'MATLAB:rmpath:DirNotFound');
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
    disp( ['    Run Path: ' params.CFGinfo.runPath])
    disp( ['    Save Path: ' params.saveDirectory])
    disp( ['    User Settings: ' params.paramsPath, '.mat' ])
    disp(' ')
    disp('--------------------------------------------------------------------------------------')

    % Locate the single cell movie
    if isempty( params.cellinfo.moviePath) 
        [cfile, cpath] = uigetfile({'*.mat';'*.tiff'});
        if cfile == 0
            error( 'singleCell: no file was selected for import'); end
        params.moviePath = [cpath, cfile];
    end

    % Import the single cell movie
    cellData = importSingleCell( params.cellinfo.moviePath);

    % Run the specific type of cell
    switch params.cellinfo.celltype
        case 'Mitosis'
            
            % initialize a mitotic cell
            myCell = MitoticCell( cellData.cell3D, params.cellinfo.lifetime, params.cellinfo.species, params.cellinfo.channelFeatures, params.cellinfo.channelsToFit, params );

            % Fit features
            myCell = myCell.fitFeatures();

        case 'Interphase'

            % Not set up yet
            error('singleCell: interphase not set up yet')

    end

end


