function paramsPath = initParams( CFG )
    % Creates the settings for running the software
    
    % set up the base configuration
    if nargin == 0; CFG='Local'; end

    % Initialize Configuration, Cell Info and Paths
    if exist( fullfile( pwd, 'initConfiguration.m') ) ~= 2
        error( ['copy initConfiguration.m to the current folder location : ', pwd] );
    else
        [params.CFG, params.CFGinfo] = initConfiguration( CFG);
    end

    % Set up the cell info for fitting
    if exist( fullfile( pwd, 'initCellInfo.m') ) ~= 2
        error( ['copy initCellInfo.m to the current folder location : ', pwd] );
    else
        params.cellinfo = initCellInfo();
    end

    % St up the paths for running single cell
    if exist( fullfile( pwd, 'initPaths.m') ) ~= 2
        error( ['copy initPaths.m to the current folder location : ', pwd] );
    else
        params = initPaths( params);
    end
    
    % Turn on diary to record all information
    c1 = onCleanup( @() eval('diary off') );
    diary( fullfile(params.saveDirectory, 'singleCell.log') );

    % Save params
    paramsPath = params.paramsPath;
    save( paramsPath, 'params', '-v7.3')
    
    % Make sure you copy the settings file from a default place into runpath

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

end
