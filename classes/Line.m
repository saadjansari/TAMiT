classdef Line < BasicElement

    properties
        startPosition
        endPosition
        length
        orientation
    end

    methods
       
        % Line {{{
        function obj = Line( startPosition, endPosition, amplitude, sigma, dim, image, props2Fit, display)
        % Line : this is the constructor function for a Line. This could be a Microtubule

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( size( image) ) || dim ~= length( startPosition) || dim ~= length(endPosition) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end
            obj = obj@BasicElement( dim, image, amplitude, sigma, props2Fit, display, 'Line');

            obj.startPosition = startPosition;
            obj.endPosition = endPosition;
            obj.length = norm( endPosition-startPosition);

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props2get)

            % sample props2get
            if nargin==1
            props2get = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            end
            
            % make sure props input matches properties defined in class
            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
            end

            % Get vector of Properties
            vec = [];
            for jProp = 1 : length( props2get)
                vec = [ vec, obj.( props2get{jProp} ) ];
            end

            % Also get a string array with property names
            vecLabels = {};
            for jProp = 1 : length( props2get)
                if numel( obj.( props2get{jProp} ) ) ~= length( obj.( props2get{jProp} ) )
                    error('getVec : property selected is a non-singleton matrix')
                end
                numRep = length( obj.( props2get{jProp} ) );
                labelRep = cell( 1, numRep);
                labelRep(:) = props2get(jProp);
                vecLabels = { vecLabels{:}, labelRep{:} };
            end

        end
        % }}}
        
        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)

            props2find = {'startPosition', 'endPosition', 'amplitude', 'sigma'};

            % find the index of start positions
            for jProp = 1 : length( props2find)
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );
                
                % Checking
                if length( idxProp) == 0
                    continue
                end
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVec: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{jProp} );
                vec( idxProp);
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}
        
        % simulateFeature {{{
        function imageOut = simulateFeature( obj, imageIn)

            if nargin < 2
                imageIn = 0*obj.image;
            end

            % Simulate a gaussian line
            if ~isempty( obj.params.idxVoxels)
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, obj.sigma, 0*imageIn, obj.params.idxVoxels.idx, obj.params.idxVoxels.X, obj.params.idxVoxels.Y, obj.params.idxVoxels.Z) );
            else
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, obj.sigma, 0*imageIn) );
            end
            imageOut = imageIn;
            obj.imageSim = imageFeat;
%             imageOut( imageFeat > imageIn) = imageFeat( imageFeat > imageIn);
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

            % Create the line to display
            line( [obj.startPosition(1) obj.endPosition(1)], [obj.startPosition(2) obj.endPosition(2)], obj.display{:} )

        end
        % }}}
        
        function obj = fillParams( obj)
            
            obj.params.idxVoxels = Line.findVoxelsNearLine( obj.startPosition, obj.endPosition, obj.image, 20);
            
        end
        
    end
    methods ( Static = true )
        
        function lineVox = findVoxelsNearLine( startPoint, endPoint, imageFind, radDilate)
            
            dim = numel( startPoint);
            if dim ~= 2 && dim ~= 3
                error('findVoxelsNearLine: dim must be 2 or 3')
            end
            
            % Draw the line in empty space
            imLine = 0*imageFind;
            len = norm( endPoint - startPoint);
            
            X = round( linspace( startPoint(1), endPoint(1), round(len) ) );
            Y = round( linspace( startPoint(2), endPoint(2), round(len) ) );
            if dim == 3
                Z = round( linspace( startPoint(3), endPoint(3), round(len) ) );
            end
            
            idxVoxels = sub2ind( size(imageFind), Y, X, Z);
            imLine( idxVoxels ) = 1;
            
            % dilate the line with a sphere of large size
            % Assume z-direction is limited (~10 pixels)
            imLine = imdilate( max( imLine, [], 3), strel('disk', radDilate) );
            imLine = repmat( imLine, 1, 1, size(imageFind, 3) );
%             imLine = imdilate( imLine, strel( 'sphere', radDilate) );
            
            % return voxel indices for the dilated line
            lineVox.idx = find( imLine(:) );
            if dim == 2
                [lineVox.Y, lineVox.X] = ind2sub( size(imLine), lineVox.idx);
            elseif dim == 3
                [lineVox.Y, lineVox.X, lineVox.Z] = ind2sub( size(imLine), lineVox.idx);
            end
            
            
        end

    end
end
