function userSettingsPath = makeSettings( CFG )
    % Creates the settings for running the software
    
    % set up the base configuration
    if nargin == 0; CFG='Local'; end
    [params.CFG, params.CFGinfo] = initConfiguration( CFG);

    % Set up the cell info for fitting
    params.cellinfo = initCellInfo();

    % Create name for save Directory. Create save directory
    [~,cellTag, ~] = fileparts( params.cellinfo.moviePath);
    resultsDirectory = [ params.CFGinfo.saveParent, filesep, 'Results' ];
    saveDirectory = [resultsDirectory, filesep, cellTag];
    warning('off', 'MATLAB:MKDIR:DirectoryExists'); addpath( resultsDirectory); mkdir( resultsDirectory, cellTag); warning('on', 'MATLAB:MKDIR:DirectoryExists')

    params.cellinfo.cellTag = cellTag;
    params.saveDirectory = saveDirectory;

    % Save this in a .mat file
    userSettingsPath = [ saveDirectory, filesep, 'settings'];
    % Make sure you copy the settings file from a default place into runpath
    save( userSettingsPath, 'params', '-v7.3')

    disp('--------------------------------------------------------------------------------------')
    disp(' ')
    disp('Paths:')
    disp(' ')
    disp( ['    Run Path: ' params.CFGinfo.runPath])
    disp( ['    Save Path: ' params.saveDirectory])
    disp( ['    User Settings: ' userSettingsPath, '.mat' ])
    disp(' ')
    disp('--------------------------------------------------------------------------------------')

end
