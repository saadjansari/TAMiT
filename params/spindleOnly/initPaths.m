function paramsCell = initPaths( params)

    numCells = length( params.cellinfo.moviePath);
    moviePaths = params.cellinfo.moviePath;
    lifetimes = params.cellinfo.lifetime;
    paramsCell = cell(1, numCells);

    for jCell = 1 : numCells

        params.cellinfo.moviePath = moviePaths{ jCell};
        params.cellinfo.lifetime = lifetimes{ jCell};

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
        params.paramsPath = [ params.saveDirectory, filesep, 'params.mat'];

        paramsCell{ jCell} = params;

    end

end
