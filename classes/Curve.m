classdef Curve < BasicElement

    properties
        startPosition
        endPosition
        length
%         orientation
    end

    methods
       
        % Line {{{
        function obj = Curve( startPosition, endPosition, amplitude, sigma, dim, props2Fit, display)
        % Line : this is the constructor function for a Line. This could be a Microtubule

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( startPosition) || dim ~= length(endPosition) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end
            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'Line');

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
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVec: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}
        
        % simulateFeature {{{
        function imageOut = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: input needed for size of image to simulate the feature in ')
            end
            imageOut = zeros( sizeImage);

            % Simulate a gaussian line
            if isfield( obj.params, 'idx') 
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, obj.sigma, imageOut, obj.params.idx, obj.params.x, obj.params.y, obj.params.z) );
            else
                imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, obj.sigma, imageOut) );
            end
            obj.imageSim = imageFeat;
            imageOut = imageFeat + imageOut;
            %imageOut( imageFeat > imageIn) = imageFeat( imageFeat > imageIn);

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                error('displayFeature: must provide axes handle to display the feature in')
            end

            % Create the line to display
            line( [obj.startPosition(1) obj.endPosition(1)], [obj.startPosition(2) obj.endPosition(2)], obj.display{:} )

        end
        % }}}
        
        % fillParams {{{
        function obj = fillParams( obj, sizeImage)
            
            obj.params.idxVoxels = Line.findVoxelsNearLine( obj.startPosition, obj.endPosition, sizeImage, 20);
            
        end
        % }}}
        
        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.props2Fit = obj.props2Fit;
            S.startPosition = obj.startPosition;
            S.endPosition = obj.endPosition;
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.display = obj.display;

        end
        % }}}
        
    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Line')
                error('incorrect type')
            end

            obj = Line( S.startPosition, S.endPosition, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);
%             obj = obj@BasicElement( S.dim, S.amplitude, S.sigma, S.props2Fit, S.display, 'Line');
%             obj.startPosition = startPosition;
%             obj.endPosition = endPosition;
%             obj.length = norm( endPosition-startPosition);

        end
        % }}}
        
        % findVoxelsNearLine {{{
        function lineVox = findVoxelsNearLine( startPoint, endPoint, sizeImage, radDilate)
            
            dim = numel( startPoint);
            if dim ~= 2 && dim ~= 3
                error('findVoxelsNearLine: dim must be 2 or 3')
            end
            
            % Draw the line in empty space
            imLine = zeros( sizeImage);
            len = norm( endPoint - startPoint);
            
            X = round( linspace( startPoint(1), endPoint(1), round(len) ) );
            Y = round( linspace( startPoint(2), endPoint(2), round(len) ) );
            if dim == 3
                Z = round( linspace( startPoint(3), endPoint(3), round(len) ) );
            end
            
            idxVoxels = sub2ind( sizeImage, Y, X, Z);
            imLine( idxVoxels ) = 1;
            
            % dilate the line with a sphere of large size
            % Assume z-direction is limited (~10 pixels)
            imLine = imdilate( max( imLine, [], 3), strel('disk', radDilate) );
            imLine = repmat( imLine, 1, 1, sizeImage(3) );
%             imLine = imdilate( imLine, strel( 'sphere', radDilate) );
            
            % return voxel indices for the dilated line
            lineVox.idx = find( imLine(:) );
            if dim == 2
                [lineVox.Y, lineVox.X] = ind2sub( sizeImage, lineVox.idx);
            elseif dim == 3
                [lineVox.Y, lineVox.X, lineVox.Z] = ind2sub( sizeImage, lineVox.idx);
            end
            
            
        end
        % }}}

    end
end
