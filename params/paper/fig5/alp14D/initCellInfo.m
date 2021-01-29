function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    type = 'Monopolar';
    species = 'Pombe';
    strain = 'alp14D';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1];

    % Lifetimes
    lifetime = { ...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 135],...
                [1 135],...
                [1 135],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 141],...
                [1 141],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 142],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 142],...
        };
 
    moviePath = {...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_001_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_001_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_001_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_001_4.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_001_5.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_002_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_002_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_002_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_003_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_003_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_003_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_004_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_004_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_004_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_005_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_005_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_006_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_006_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_006_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_007_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_007_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_007_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_008_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_008_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_008_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_009_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_009_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_009_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_010_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_010_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_010_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_011_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_011_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_011_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_011_4.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_012_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_012_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_012_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_012_4.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_012_5.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_013_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_013_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_013_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_014_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_014_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_014_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_015_1.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_015_2.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_015_3.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_015_4.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_015_5.mat',...
         'Paper/Monopolar/16-10-26 (Alp14D dynamics random with Piezo)/863_150ms_016_1.mat',...
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