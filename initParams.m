function paramsPath = initParams( opts)
    % Creates the settings for running the software
    
    % ---------------------  PATHS SPEC  -------------------------
    % Run path: Location of main.m
    % Save path: Path to results folder for storage
    % Local
    exec_loc.local.runpath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';
    exec_loc.local.savepath = fullfile( exec_loc.local.runpath, 'Results');
    % Summit
    exec_loc.summit.runpath = '/projects/saan8193/ImageAnalysis/SingleCell';
    exec_loc.summit.savepath = '/scratch/summit/saan8193/SingleCell';
    % ------------------------------------------------------------

    
    % Initialize Configuration, Cell Info and Paths
    if exist( fullfile( pwd, 'initConfiguration.m') ) ~= 2
        error( ['copy initConfiguration.m to the current folder location : ', pwd] );
    else
        params = initConfiguration( opts);
    end

    % Set up the cell info for fitting
    if exist( fullfile( pwd, 'initCellInfo.m') ) ~= 2
        error( ['copy initCellInfo.m to the current folder location : ', pwd] );
    else
        params.cellInfo = initCellInfo( opts);
    end

    % Set up the paths for running single cell
    params = initPaths( params, exec_loc);

    numCells = length( params);
    paramsRef = params;

    for jCell = 1 : numCells

        params = paramsRef{ jCell};

        % Save params
        paramsPath{ jCell} = params.paramsPath;
        save( paramsPath{ jCell}, 'params', '-v7.3')

    end

    function paramsCell = initPaths( params, exec_loc)

        numCells = params.cellInfo.numCells;
        moviePaths = params.cellInfo.moviePath;
        lifetimes = params.cellInfo.lifetime;
        paramsCell = cell(1, numCells);

        % Assign paths
        % runPath : location where the main.m file will be run from.
        % saveParent : the parent directory where the results folder will reside
        switch params.LOC
            case 'Local'
                runPath = exec_loc.local.runpath;
                saveParent =  exec_loc.local.savepath;
            case 'Summit'
                runPath = exec_loc.summit.runpath;
                saveParent = exec_loc.summit.savepath;
        end
        
        params.runPath = runPath;
        params.saveParent = saveParent;
        if exist( params.saveParent) ~=7
            mkdir( params.saveParent)
        end

        % Get string for current datetime in YYMMDD_HHMM format
        timeTag = datestr( now, 'yymmdd_HHMM');

        for jCell = 1 : numCells

            params.cellInfo.moviePath = moviePaths{ jCell};
            params.cellInfo.lifetime = lifetimes{ jCell};

            % Get Cell Tag
            [~,cellTag, ~] = fileparts( params.cellInfo.moviePath);

            % Create string for directory name where data from this particular
            saveDirectory = [ params.saveParent, filesep, timeTag, '__', cellTag];
            mkdir( saveDirectory);
            addpath( saveDirectory); 

            % Save directory path in params structure
            params.cellInfo.cellTag = cellTag;
            params.saveDirectory = saveDirectory;
            params.paramsPath = [ params.saveDirectory, filesep, 'params.mat'];
            paramsCell{ jCell} = params;

        end

    end


end
