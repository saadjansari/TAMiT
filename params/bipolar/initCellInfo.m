function cellInfo = initCellInfo( opts)
    % Specifies information about the cells of interest

    % Path to parent folder containing segmented cell movies
    movieParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets';
 
    % Relative path to segment cell movies
    moviePath = {...
         'LukeTest/1149_100R_100G_25_deg_001_A.mat',...
%          'LukeTest/1149_100R_100G_25_deg_001_B.mat',...
%          'LukeTest/1149_100R_100G_25_deg_001_C.mat',...
%          'LukeTest/1149_100R_100G_25_deg_001_D.mat',...
%          'LukeTest/1149_100R_100G_25_deg_001_E.mat',...
    };

    % Frames to analyze for each cell
    lifetime = { ...
                [100 102],...
%                 [],...
%                 [],...
%                 [],...
%                 [],...
        };

    % Cell Type : 'Mitosis' / 'Monopolar'
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
    species = 'Pombe';
    strain = 'WT_sid4';

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
