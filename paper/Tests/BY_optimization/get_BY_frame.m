function imData = get_BY_frame( cellpath, cellname, time)
    % Load BY movies and return some number of frames for each cell
    addpath('/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/functions/external/bfmatlab')
        
    reader = bfGetReader( strjoin({cellpath, cellname }, filesep) );
    metaData = reader.getMetadataStore();
    Meta.numVoxelsX = metaData.getPixelsSizeX(0).getValue();
    Meta.numVoxelsY = metaData.getPixelsSizeY(0).getValue();
    Meta.numVoxelsZ = metaData.getPixelsSizeZ(0).getValue();
    Meta.numTimes = metaData.getPixelsSizeT(0).getValue();
    Meta.numChannels = metaData.getPixelsSizeC(0).getValue();

    % get physical sizes of voxels
    Meta.sizeVoxelsX = double( metaData.getPixelsPhysicalSizeX(0).value() );
    Meta.sizeVoxelsY = double( metaData.getPixelsPhysicalSizeY(0).value() );
    Meta.sizeVoxelsZ = double( metaData.getPixelsPhysicalSizeZ(0).value() );

    % Get random time indices
    fprintf('Cell = %s \nTimes =\n', cellname )
    disp(times)

    % Get data type of the images
    im_temp = bfGetPlane( reader, reader.getIndex( 0, 0, 0) + 1);

    % Pre-allocate 3D array for storing the image planes:
    imData = zeros( Meta.numVoxelsY, Meta.numVoxelsX, Meta.numVoxelsZ, class(im_temp) );

    % Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
    % All indices are expected to be 1-based
    reader.setSeries(0);
    for jZ = 1 : Meta.numVoxelsZ
        iPlane = reader.getIndex(jZ-1, 0, time-1) + 1;
        imData(:,:,jZ) = bfGetPlane(reader, iPlane);
    end
end