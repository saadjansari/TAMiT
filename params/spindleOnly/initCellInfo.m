function params = initCellInfo()
    % Initializes the cell info

    % Cell Type : 'Mitosis' or 'Interphase'
    celltype = 'Mitosis';

    % Cell Species
    species = 'Pombe';

    % Cell Strain
    strain = 'gtb1-K5A';

    % Lifetimes
    lifetime = { ...
                [35 36], ...
                [31 32], ...
                [10 15], ...
        };

    % Channel Features ( what features do you want to fit?)
    channelFeatures = {'Microtubule', 'Cut7'};

    % Channels to Fit ( [] = all ) (length must match the length of channelFeatures) (order must match the order of features in channelFeatures)
    channelsToFit = [2, 1];
 
    % Cell Movie Location
    movieParentLocal = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets';
    movieParentSummit = '/projects/saan8193/ImageAnalysis/FY Datasets';
    movieParentRumor = '/projects/saan8193/ImageAnalysis/FY Datasets';
    moviePath = {...
        '1095_50msG_100msR_7Z_005_cells/1095_50msG_100msR_7Z_005_7.mat', ...
        '1095_50msG_100msR_7Z_005_cells/1095_50msG_100msR_7Z_005_26.mat', ...
        '1095_50msG_100msR_7Z_005_cells/1095_50msG_100msR_7Z_005_33.mat', ...
    };
    moviePathFull = cellfun( @(x) fullfile( movieParentLocal, x), moviePath, 'UniformOutput', false); 
    
    params.celltype = celltype;
    params.species = species;
    params.strain = strain;
    params.lifetime = lifetime;
    params.channelFeatures = channelFeatures;
    params.channelsToFit = channelsToFit;
    params.moviePath = moviePathFull;

end
