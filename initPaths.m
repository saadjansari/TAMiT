function paramsCell = initPaths( params)

    numCells = params.cellInfo.numCells;
    moviePaths = params.cellInfo.moviePath;
    lifetimes = params.cellInfo.lifetime;
    paramsCell = cell(1, numCells);

    % Assign paths
    % runPath : location where the main.m file will be run from.
    % saveParent : the parent directory where the results folder will reside
    switch params.LOC
        case 'Local'
            runPath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';
            saveParent = [runPath,filesep,'Results'];

        case 'Summit'
            runPath = '/projects/saan8193/ImageAnalysis/SingleCell';
            saveParent = '/scratch/summit/saan8193/SingleCell';
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
