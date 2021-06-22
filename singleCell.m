function features = singleCell( paramsPath)
    % Input: paramsPath is a path to a .mat file

    clc; close all;
    clearvars -except paramsPath 

    % Make sure we have the correct paths
    addpath( pwd);
    addpath( genpath( [pwd, filesep, 'functions'] ) );
    addpath( [pwd, filesep, 'classes'] );
    addpath( [pwd, filesep, 'classes/accessory'] );

    % Load the params
    load( paramsPath); 

    disp('--------------------------------------------------------------------------------------')
    disp( ['Configuration: ' params.CFG])
    disp(' ')
    disp('Paths:')
    disp(' ')
    disp( ['    Run Path: ' params.runPath])
    disp( ['    Save Path: ' params.saveDirectory])
    disp( ['    User Settings: ' params.paramsPath])
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

    % Start parallel pool
    if params.fit.useParallel
        pc = parcluster('local');
        if strcmp(params.LOC, 'Summit')
            % explicitly set the JobStorageLocation to a temp directory
            pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
            parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))
        elseif strcmp(params.LOC, 'Local')
            parpool(pc);
        else
            error('boop')
        end
    end

    % Run the specific type of cell
    switch params.cellInfo.type
        case 'Mitosis'
            fname = 'MitoticCell';
        case 'Monopolar'
            fname = 'MonopolarCell';
        case 'Interphase'
            fname = 'InterphaseCell';
        case 'MitosisBud'
            fname = 'MitoticCellBud';
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


