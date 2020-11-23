function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'Monopolar';
    species = 'Pombe';
    strain = 'WT';

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
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 141],...
                [1 141],...
                [1 141],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 140],...
                [1 140],...
                [1 70],...
                [1 141],...
                [1 80],...
                [1 80],...
                [1 80],...
                [1 80],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 139],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 139],...
                [1 139],...
                [1 140],...
                [1 140],...
                [1 140],...
                [1 139],...
                [1 139],...
        };
 
    moviePath = {...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_001_6.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_002_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_002_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_002_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_002_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_002_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_003_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_003_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_003_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_004_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_004_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_004_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_005_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_005_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_006_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_006_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_006_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_006_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_006_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_007_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_007_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_007_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_008_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_008_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_008_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_009_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_009_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_009_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_009_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_009_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_010_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_010_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_010_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_010_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_010_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_011_6.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_012_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_012_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_012_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_012_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_013_6.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_014_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_014_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_014_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_014_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_014_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_015_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_015_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_015_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_015_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_015_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_016_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_016_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_016_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_016_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_016_5.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_017_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_017_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_018_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_019_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_020_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_020_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_020_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_020_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_021_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_021_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_021_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_021_4.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_022_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_022_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_023_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_023_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_024_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_024_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_024_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_025_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_025_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_026_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_026_2.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_026_3.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_027_1.mat',...
         'Paper/Monopolar/16-10-20 (WT dynamics random with Piezo)/850_150ms_027_2.mat',...
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
