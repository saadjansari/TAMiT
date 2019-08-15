function features = singleCell( settingsPath)
    % Input: settingsPath is a path to a .mat file

    % Load the settings
    settings = load( settingsPath); settings = settings.params;

    % Locate the single cell movie
    if isempty( settings.cellinfo.moviePath) 
        [cfile, cpath] = uigetfile({'*.mat';'*.tiff'});
        if cfile == 0
            error( 'singleCell: no file was selected for import'); end
        settings.moviePath = [cpath, cfile];
    end

    % Import the single cell movie
    cellData = importSingleCell( settings.cellinfo.moviePath);

    % Run the specific type of cell
    switch settings.cellinfo.celltype
        case 'Mitosis'
            
%             parpool;

            % initialize a mitotic cell
            myCell = MitoticCell( cellData.cell3D, settings.cellinfo.lifetime, settings.cellinfo.species, settings.cellinfo.channelFeatures, settings );

            % Play the movie
%             myCell.playMovie( 1 );

            % Find features
            myCell = myCell.fitFeatures();

%             delete(gcp);

        case 'Interphase'

            % Not set up yet
            error('singleCell: interphase not set up yet')

    end

end


