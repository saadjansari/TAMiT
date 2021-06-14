function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'Mitosis';
    species = 'Pombe';
    strain = 'WT_sid4';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Sid4', 'Cut7'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1,2];

    % Lifetimes
    lifetime = { ...
                [],...
        };
 
    moviePath = {...
         'sid1_cut7/1032_100msR_50msG_7Z_106_cells/1032_100msR_50msG_7Z_106_1.mat',...
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
