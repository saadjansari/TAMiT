function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'MitosisBud';
    species = 'BuddingYeast';
    strain = 'test';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1];

    % Lifetimes
    lifetime = { ...
                [], ...
                [], ...
                [], ...
                [], ...
                [], ...
                [], ...
        };
 
    moviePath = {...
        '1.10.18 37C 1.tif', ...
        '1.10.18 37C 2.tif', ...
        '1.10.18 37C 4 curved MTs.tif', ...
        '1.24.18 25C 3.tif', ...
        '1.24.18 37C 6.tif', ...
        '2.7.18 37C 2,3 curved MT.tif', ...
    };

    % Cell Movie Location
    movieParentLocal = '/Users/saadjansari/Documents/Projects/ImageAnalysis/BY Datasets';
    movieParentSummit = '/projects/saan8193/ImageAnalysis/BY Datasets';
    movieParentRumor = '/projects/saan8193/ImageAnalysis/BY Datasets';
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
