function cellInfo = initCellInfo( opts)
    % Specifies information about the cells of interest

    % Path to parent folder containing segmented cell movies
    movieParent = [pwd, filesep, 'demos', filesep, 'FissionYeastBipolar'];
 
    % Relative path to segment cell movies
    moviePath = {...
          'cell_data.mat',...
    };

    % Frames to analyze for each cell
    lifetime = { ...
                [80 85],...
        };

    % Cell Type : 'Mitosis' / 'Monopolar' / 'MitosisBud'
    type = 'Mitosis';

    % Channels to Fit (list of channel numbers)
    % Example: [1], [1,2], [] (empty brackets denote all channels)
    % Note: length must match the length of channelFeatures. Order must also match the order of features in channelFeatures)
    channelsToFit = [1];

    % Features to fit
    % Current available features:
    % 'Microtubule', 'Sid4', 'Cut7', 'KC'
    channelFeatures = {'Microtubule'};

    % Additional book-keeping information
    species = 'FY';
    strain = 'Demo';

    % Cell Movie Location
    movieParentLocal = movieParent;
    movieParentSummit = '/projects/saan8193/ImageAnalysis/FY Datasets';
    movieParentRumor = '/projects/saan8193/ImageAnalysis/FY Datasets';
    switch opts.LOC
        case 'Local'
            movieParent = movieParentLocal;
        case 'Summit'
            movieParent = movieParentSummit;
        case 'Rumor'
            movieParent = movieParentRumor;
    end
    
    numCells = length( moviePath);
    % Construct full movie paths
    for jCell = 1 : numCells
        moviePathFull{jCell} = fullfile( movieParent, moviePath{jCell});
    end

    
    cellInfo.type = type;
    cellInfo.species = species;
    cellInfo.strain = strain;
    cellInfo.lifetime = lifetime;
    cellInfo.channelFeatures = channelFeatures;
    cellInfo.channelsToFit = channelsToFit;
    cellInfo.moviePath = moviePathFull;
    cellInfo.numCells = numCells;

end
