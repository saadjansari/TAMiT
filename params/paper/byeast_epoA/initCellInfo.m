function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'MitosisBud';
    species = 'BuddingYeast';
    strain = 'EpoA';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1];

    % Lifetimes
    lifetime = { ...
                %[], ...
                [], ...
                %[], ...
                %[], ...
                %[], ...
                %[], ...
        };
 
    moviePath = {...
        %'EpoA MTs/1.25.17 2306 EpoA002-1.tif', ...
        'EpoA MTs/1.25.17 2306 EpoA003-1.tif', ...
        %'EpoA MTs/1.25.17 2306 EpoA004-1.tif', ...
        %'EpoA MTs/1.25.17 2306 EpoA006-1.tif', ...
        %'EpoA MTs/2.16.17 2307 10um EpoA003-1.tif', ...
        %'EpoA MTs/2.16.17 2307 10um EpoA004-1.tif', ...
        %'EpoA MTs/6.2.17 2306 10uM EpoA002-1.tif', ...
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
