function userSettingsPath = makeSettings()

    % Cell Movie Location
    moviePath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/Test Data Sets/MT interphase to monopolar to bipolar/857_100msR_7Z_001_cells/857_100msR_7Z_001_28.mat';
    % Cell Type : 'Mitosis' or 'Interphase'
    celltype = 'Mitosis';
    % Cell Species
    species = 'Pombe';
    % Cell Strain
    strain = 'wild-type';
    % Times to Analyze: must be within the time range of the cell movie
    lifetime = [208 210]; % [150 155]

    % Channel Features
%     channelFeatures = {'Microtubule', 'Kinetochore'};
    channelFeatures = {'Microtubule'};

    % Server Info : 'Local' or 'Summit'
    server = 'Local';

    % Paths
    % runPath : location where the main.m file will be run from.
    % saveParent : the parent directory where the results folder will reside
    if strcmp( server, 'Summit') 

        runPath = '/projects/saan8193/ImageAnalysis/SingleCell';
        saveParent= '/scratch/summit/saan8193/ImageAnalysis/SingleCell'

    elseif strcmp( server, 'Local')

        runPath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';
        saveParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';

    end

    % Save Directory
    [~,cellTag, ~] = fileparts( moviePath);
    resultsDirectory = [ saveParent, filesep, 'Results' ];
    saveDirectory = [resultsDirectory, filesep, cellTag];

    warning('off', 'MATLAB:MKDIR:DirectoryExists')
    mkdir( resultsDirectory, cellTag);
    warning('on', 'MATLAB:MKDIR:DirectoryExists')

%     interactive = 1; % 1 for real-time plots, 0 to save the data for plotting

    % Save this in a .mat file
    userSettingsPath = [ saveDirectory, filesep, 'settings'];
    % Make sure you copy the settings file from a default place into runpath
    save( userSettingsPath, '-v7.3')

    disp('--------------------------------------------------------------------------------------')
    disp(' ')
    disp('Paths:')
    disp(' ')
    disp( ['    Run Path: ' runPath])
    disp( ['    Save Path: ' saveDirectory])
    disp( ['    User Settings: ' userSettingsPath, '.mat' ])
    disp(' ')
    disp('--------------------------------------------------------------------------------------')

end
