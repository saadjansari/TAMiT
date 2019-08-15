function params = initPaths( params)

    % Get Cell Tag
    [~,cellTag, ~] = fileparts( params.cellinfo.moviePath);
    
    % Get string for current datetime in YYMMDD_HHMM format
    timeTag = datestr( now, 'yymmdd_HHMM');
    
    % Create string for directory name where data from this particular
    % analysis will be saved
    saveDirectory = [ params.CFGinfo.saveParent, filesep, timeTag, '__', cellTag];

    % Create save directory and add it to the path
    if exist( params.CFGinfo.saveParent) ~=7
        mkdir( params.CFGinfo.saveParent)
    end
    mkdir( saveDirectory);
    addpath( saveDirectory); 
    
    % Save directory path in params structure
    params.cellinfo.cellTag = cellTag;
    params.saveDirectory = saveDirectory;
    
    params.paramsPath = [ params.saveDirectory, filesep, 'params'];

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