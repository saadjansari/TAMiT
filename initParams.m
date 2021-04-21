function paramsPath = initParams( opts)
    % Creates the settings for running the software
    %
    
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
    if exist( fullfile( pwd, 'initPaths.m') ) ~= 2
        error( ['copy initPaths.m to the current folder location : ', pwd] );
    else
        params = initPaths( params);
    end

    numCells = length( params);
    paramsRef = params;

    for jCell = 1 : numCells

        params = paramsRef{ jCell};

        % Save params
        paramsPath{ jCell} = params.paramsPath;
        save( paramsPath{ jCell}, 'params', '-v7.3')

    end


end
