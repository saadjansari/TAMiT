function params = initCellInfo()
    % Initializes the cell info
    

    % Cell Type : 'Mitosis' or 'Interphase'
    celltype = 'Mitosis';

    % Cell Species
    species = 'Pombe';

    % Cell Strain
    
    strain = 'wild-type';

    % Times to Analyze: must be within the time range of the cell movie
    lifetime = [150 155];
    %lifetime = [208 210]; 

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule', 'Kinetochore'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [1 2];

    % Cell Movie Location
    moviePath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/998_150msR_100G_trig_7Z_001_cells/998_150msR_100G_trig_7Z_001_24.mat'; 
%     moviePath = '/projects/saan8193/ImageAnalysis/FY Datasets/998_150msR_100G_trig_7Z_001_cells/998_150msR_100G_trig_7Z_001_24.mat';     
    
    params.celltype = celltype;
    params.species = species;
    params.strain = strain;
    params.lifetime = lifetime;
    params.channelFeatures = channelFeatures;
    params.channelsToFit = channelsToFit;
    params.moviePath = moviePath;

end
