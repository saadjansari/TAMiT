function features = singleCell( settingsPath)
    % Input: settingsPath is a path to a .mat file

    % Load the settings
    settings = load( settingsPath);

    % Locate the single cell movie
    if isempty( settings.moviePath) 
        [cfile, cpath] = uigetfile({'*.mat';'*.tiff'});
        if cfile == 0
            error( 'singleCell: no file was selected for import'); end
        settings.moviePath = [cpath, cfile];
    end

    % Import the single cell movie
    cellData = importSingleCell( settings.moviePath);

    % Run the specific type of cell
    switch settings.celltype
        case 'Mitosis'
            
            % initialize a mitotic cell
            myCell = MitoticCell( cellData.cell3D, settings.lifetime, settings.species, settings.channelFeatures, settings );

            % Play the movie
%             myCell.playMovie( 1 );

            % Find features
            myCell = myCell.fitFeatures();

        case 'Interphase'

            % Not set up yet
            error('singleCell: interphase not set up yet')

    end

end


