classdef Spot < BasicElement

    properties
        position
    end

    methods ( Access = public )
       
        % Spot {{{
        function obj = Spot( position, amplitude, sigma, dim, props2Fit, display)
        % Spot : this is the constructor function for a Spot. This will be used for SPBs, MTOCs, and Kinetochores

            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'Spot');
            obj.position = position;

        end
        % }}}

        % simulateFeature {{{
        function imageOut = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: must provide size of box to simulate feature in')
            end

            % Simulate a gaussian spot
            imageOut = zeros( sizeImage);
            if isfield( 'idxVoxels', obj.params)
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianPoint3D( obj.position, obj.sigma, imageOut, obj.params.idxVoxels.idx, obj.params.idxVoxels.X, obj.params.idxVoxels.Y, obj.params.idxVoxels.Z) );
            else
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianPoint3D( obj.position, obj.sigma, imageOut) );
            end
            obj.imageSim = imageFeat;
%             imageOut( imageOut < imageFeat) = imageFeat( imageOut < imageFeat);
            imageOut = imageFeat + imageOut;

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.image, [], 3) ); colormap gray; axis equal;
            end

            plot( obj.position(1), obj.position(2), obj.display{:} ) 

        end
        % }}}

        % BleachSpot {{{
        function imgBleached = BleachSpot( obj, imgIn, bleachRadius)
            % Bleaches the Spot location with a circle of zeros of radius bleachRad
            
            if length( size( imgIn) ) ~= 2
                error('BleachIt : dimensionality of the input image must be 2'); end

            x0 = round( obj.position(1) ); y0 = round( obj.position(2) );
            
            imgBleached = logical( 0* imgIn );
            imgBleached( y0, x0) = 1;
            imgBleached = imdilate( imgBleached, strel('disk', bleachRadius) ); % circle with ones
            
            imgBleached = double( ~imgBleached); % circle with zeros
            
            imgBleached = imgBleached .* imgIn;

        end
        % }}}

        % fillparams {{{
        function obj = fillParams( obj)
            
            obj.params.idxVoxels = Spot.findVoxelsNearSpot( obj.position, obj.image, 15);
            
        end
        % }}}
        
        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.props2Fit = obj.props2Fit;
            S.position = obj.position;
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.display = obj.display;

        end
        % }}}
        
    end
    
    methods (Static = true)

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Spot')
                error('incorrect type')
            end

            obj = Spot( S.position, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);
%             obj = obj@BasicElement( S.dim, S.amplitude, S.sigma, S.props2Fit, S.display, 'Spot');

        end
        % }}}
        
        % findVoxelsNearSpot {{{
        function spotVox = findVoxelsNearSpot( pos, imageFind, radDilate)
            
            dim = numel( pos );
            if dim ~= 2 && dim ~= 3
                error('findVoxelsNearSpot: dim must be 2 or 3')
            end
            
            % Draw the spot in empty space
            imSpot = 0*imageFind;
            if dim == 2
                imSpot( round(pos(2)), round( pos(1) ) ) = 1;
            elseif dim == 3
                imSpot( round(pos(2)), round( pos(1) ), round( pos(3) ) ) = 1;
            end
            
            % dilate the line with a sphere of large size
%             imSpot = imdilate( imSpot, strel( 'sphere', radDilate) );
            imSpot = imdilate( max( imSpot, [], 3), strel('disk', radDilate) );
            imSpot = repmat( imSpot, 1, 1, size(imageFind, 3) );
            
            % return voxel indices for the dilated line
            spotVox.idx = find( imSpot(:) );
            if dim == 2
                [spotVox.Y, spotVox.X] = ind2sub( size(imSpot), spotVox.idx);
            elseif dim == 3
                [spotVox.Y, spotVox.X, spotVox.Z] = ind2sub( size(imSpot), spotVox.idx);
            end
            
            
        end
        % }}}
        
    end

end