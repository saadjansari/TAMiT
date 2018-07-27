    % This block of code takes the user selected ND2 file and imports it,
    % extracts metadata
               
    % Prompt user to select the appropriate file.
    [ filename, pathname] = uigetfile( '*.*', 'Select a file');      
    if filename == 0
        error('File must be selected for import');
    end
    [~, onlyName, ext] = fileparts(filename);
    if ~strcmp( ext,'.nd2')
        error('Must select a .nd2 file')
    end

    % Loop through files and extract metadata and image data
    reader = bfGetReader([pathname, filename]); % Creates a Bioformats reader objects
    metaData = reader.getMetadataStore();   % Use class function to acess metadata of image
    
    % Extract key parameters from  metadata
    % get number of voxels in all the dimensions
    Meta.numVoxelsX = metaData.getPixelsSizeX(0).getValue();
    Meta.numVoxelsY = metaData.getPixelsSizeY(0).getValue();
    Meta.numVoxelsZ = metaData.getPixelsSizeZ(0).getValue();
    Meta.numTimes = metaData.getPixelsSizeT(0).getValue();
    Meta.numChannels = metaData.getPixelsSizeC(0).getValue();
    
    % get physical sizes of voxels
    Meta.sizeVoxelsX = metaData.getPixelsPhysicalSizeX(0).value();
    Meta.sizeVoxelsY = metaData.getPixelsPhysicalSizeY(0).value();
    Meta.sizeVoxelsZ = metaData.getPixelsPhysicalSizeZ(0).value();
    
    % Get other useful information
    Meta.numPlanes = metaData.getPlaneCount(0);
    Meta.wavelength = double(metaData.getChannelEmissionWavelength(0,0).value);
    Meta.exposureTime = double(metaData.getPlaneExposureTime(0,0).value);
    Meta.fileName = onlyName;
    Meta.fileExtension = ext; 
    
    % Pre-allocate 5D array for storing the image planes:
    % array dimensions are arranged as ( x, y, z, t, c)
    imData = zeros( Meta.numVoxelsX, Meta.numVoxelsY, Meta.numVoxelsZ, ...
        Meta.numTimes, Meta.numChannels,'uint16');
    
    PlaneTimes = zeros( Meta.numVoxelsZ, Meta.numTimes, Meta.numChannels);
    
    % Add planes to imData for each channel, for every frame, for every z
    % plane.
    for jChannel = 1 : Meta.numChannels
        
        for jTime = 1 : Meta.numTimes
            
            for jZ = 1 : Meta.numVoxelsZ
                
                % get index of plane with specific z, t, and c
                jPlane = reader.getIndex( jZ - 1, jChannel - 1, jTime - 1) + 1;
                
                % get image plane corresponding to the index above
                imData( :, :, jZ, jTime, jChannel) = bfGetPlane(reader, jPlane);
        
                % get time correponding to this plane
                
                % THERE IS AN ISSUE HERE WITH RECORDING PLANE TIMES!!!

                try
                PlaneTimes( jZ, jTime, jChannel) = metaData.getPlaneDeltaT( 0, jPlane).value; 
                end
                
            end
            
        end
        
    end