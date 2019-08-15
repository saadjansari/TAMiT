function paramsPath = initParams( CFG )
    % Creates the settings for running the software
    
    % set up the base configuration
    if nargin == 0; CFG='Local'; end
    [params.CFG, params.CFGinfo] = initConfiguration( CFG);

    % Set up the cell info for fitting
    params.cellinfo = initCellInfo();

    % St up the paths for running single cell
    params = initPaths( params);
    
    % Save params
    paramsPath = params.paramsPath;
    save( paramsPath, 'params', '-v7.3')
    
    % Make sure you copy the settings file from a default place into runpath

end