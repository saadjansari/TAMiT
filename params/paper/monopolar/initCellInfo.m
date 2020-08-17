function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'Monopolar';
    species = 'Pombe';
    strain = 'cut11-7';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1];

    % Lifetimes
    lifetime = { ...
                [],...
        };
 
    moviePath = {...
         'Test/Monopolar/cut7_24__170601_klp6d/887_150msR_100msG_trig_7Z_cells/887_150msR_100msG_trig_7Z_001_1.mat', ...
    };

    % Cell Movie Location
    movieParentLocal = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets';
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
