function myCell = TAMiT_cell( paramsPath)
    % Input: paramsPath is a path to a .mat file

    close all;
    clearvars -except paramsPath 

    % Make sure we have the correct paths
    addpath( pwd);
%     addpath( genpath( [pwd, filesep, 'functions'] ) );
    addpath( [pwd, filesep, 'src'] );
    addpath( genpath([pwd, filesep, 'src/accessory']) );
    addpath( genpath([pwd, filesep, 'src/cells']) );
    addpath( genpath([pwd, filesep, 'src/features']) );
    
    % Load the params
    load( paramsPath); 
 
    % Print preamble
    disp('--------------------------------------------------------')
    disp('--------------------------------------------------------')
    disp('    _______       __  __ _ _______ ')
    disp('   |__   __|/\   |  \/  (_)__   __|')
    disp('      | |  /  \  | \  / |_   | |   ')
    disp('      | | / /\ \ | |\/| | |  | |   ')
    disp('      | |/ ____ \| |  | | |  | |   ')
    disp('      |_/_/    \_\_|  |_|_|  |_|   ')
    disp(' ')
    disp('Toolkit for Automated Microtubule Tracking')
    fprintf('\n\n')
    disp(fileread('LICENSE'))
    disp('--------------------------------------------------------')
    disp('--------------- C O N F I G U R A T I O N --------------')
    disp('--------------------------------------------------------')
    disp(' ')
    disp( ['Display Mode: ' params.DisplayState])
    disp( ['Run Path: ' params.runPath])
    disp( ['Save Path: ' params.saveDirectory])
    disp( ['Params File Path: ' params.paramsPath])
    disp(' ')
    
    % Locate the single cell movie
    if isempty( params.cellInfo.moviePath) 
        [cfile, cpath] = uigetfile({'*.mat';'*.tiff'});
        if cfile == 0
            error( 'TAMiT_cell: no file was selected for import'); end
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
    
    fprintf('--------------------------------------------------------\n')
    fprintf('-------------------  S U C C E S S !!!  ----------------\n')
    fprintf('--------------------------------------------------------\n\n\n')


end


