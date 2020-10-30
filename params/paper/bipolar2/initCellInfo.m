function cellInfo = initCellInfo( opts)
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase' or 'Monopolar'
    %type = 'Interphase';
    type = 'Mitosis';
    species = 'Pombe';
    strain = 'cut11-7';

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1];

    % Lifetimes
    lifetime = { ...
        [41 173],  
        [67 173],  
        [10 115],  
        [22 148],  
        [1 63],  
        [1 70],  
        [1 88],  
        [1 97],  
        [11 137],  
        [17 159],  
        [3 133],
        [19 159],  
        [1 130],  
        [15 114],  
        [13 161],  
        [1 113],  
        [1 163],  
        [112 239],  
        [12 162],  
        [61 239],  
        [1 127],  
        [1 101],  
        };
 
    moviePath = {...
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_001_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_001_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_002_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_002_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_003_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_003_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_003_3.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_003_4.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_004_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_005_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_005_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_005_3.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_005_4.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_006_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_007_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_007_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_009_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_010_1.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_010_2.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_010_3.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_010_4.mat',
         'Paper/Bipolar/wildtype_180201and2/979_150msR_100G_trig_7Z_010_5.mat',
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
