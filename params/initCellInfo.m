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

    % Channel Features (array size must match dim5 of cell movie)
    channelFeatures = {'Microtubule', 'Kinetochore'};

    % Cell Movie Location
    moviePath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/998_150msR_100G_trig_7Z_001_cells/998_150msR_100G_trig_7Z_001_24.mat'; 

    %moviePath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/Test Data Sets/MT interphase to monopolar to bipolar/857_100msR_7Z_001_cells/857_100msR_7Z_001_28.mat';
    
    params.celltype = celltype;
    params.species = species;
    params.strain = strain;
    params.lifetime = lifetime;
    params.channelFeatures = channelFeatures;
    params.moviePath = moviePath;

end
