classdef Line < BasicElement

    properties
        startPosition
        endPosition
        length
%         orientation
    end

    methods
       
        % Line {{{
        function obj = Line( startPosition, endPosition, amplitude, sigma, dim, props2Fit, display)
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
                if isempty( idxProp)
                    continue
                end
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
                if obj.dim == 2
                    imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine2D( obj.startPosition, obj.endPosition, ...
                        obj.sigma, imageOut, obj.params.idx, obj.params.x, obj.params.y) );
                elseif obj.dim == 3
                    imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, ...
                        obj.sigma, imageOut, obj.params.idx, obj.params.x, obj.params.y, obj.params.z) );
                end
            else
                if obj.dim == 2
                    imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine2D( obj.startPosition, obj.endPosition, obj.sigma, imageOut) );
                elseif obj.dim == 3
                    imageFeat = obj.amplitude * mat2gray( Cell.drawGaussianLine3D( obj.startPosition, obj.endPosition, obj.sigma, imageOut) );
                end
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
        
        % displayFeatureXZ {{{
        function ax = displayFeatureXZ( obj, ax)

            if nargin < 2
                error('displayFeatureXZ: must provide axes handle to display the feature in')
            end
            
            if obj.dim == 2
                error('displayFeatureXZ: must be 3-dimensional')
            end

            % Create the line to display
            line( [obj.startPosition(3) obj.endPosition(3)], [obj.startPosition(1) obj.endPosition(1)], obj.display{:} )
        end
        % }}}
        
        % displayFeatureXZ {{{
        function ax = displayFeatureYZ( obj, ax)

            if nargin < 2
                error('displayFeatureYZ: must provide axes handle to display the feature in')
            end
            
            if obj.dim == 2
                error('displayFeatureYZ: must be 3-dimensional')
            end

            % Create the line to display
            line( [obj.startPosition(3) obj.endPosition(3)], [obj.startPosition(2) obj.endPosition(2)], obj.display{:} )

        end
        % }}}
        
        % fillParams {{{
        function obj = fillParams( obj, sizeImage)
            
            obj.params = Line.findVoxelsNearLine( obj.startPosition, obj.endPosition, sizeImage, 10);
%             [obj.params.y, obj.params.x, obj.params.z] = ind2sub( sizeImage, obj.params.idx);
            
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
        
        % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
            % Get 2D projection of feature
            
            % Check object dimensionality
            if obj.dim == 2
                warning('object dimensionality is already 2')
            end
            
            obj.startPosition = obj.startPosition(1:2);
            obj.endPosition = obj.endPosition(1:2);
            
        end
        % }}}
        
        % GetProjection3DSpecific {{{
        function obj = GetProjection3DSpecific( obj)
            % Get 3D projection of feature
            
            % Check object dimensionality
            if obj.dim == 3
                warning('object dimensionality is already 3')
            end
            
            obj.startPosition(3) = 1;
            obj.endPosition(3) = 1;
            
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            
            obj.startPosition(1:2) = obj2D.startPosition(1:2);
            obj.endPosition(1:2) = obj2D.endPosition(1:2);
            obj.sigma(1:2) = obj2D.sigma(1:2);
            obj.amplitude = obj2D.amplitude;
            
        end
        % }}}
        
        % GetLength {{{
        function len = GetLength( obj)
            
            len = norm( obj.startPosition - obj.endPosition);
            obj.length = len;
            
        end
        % }}}
        
        % GetOrientation {{{
        function orientationXY = GetOrientation( obj)

            orientationXY = atan( obj.endPosition(2)-obj.startPosition(2) ./ obj.endPosition(1)-obj.startPosition(1) );

        end
        % }}}
        
        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
            % Reduce the length of the line until it lies inside the mask
            
            % Simulate image
            imFeat = obj.simulateFeature( size(mask) );
            
            % if escapes from mask, shorten until inside mask
            outside = Methods.CheckEscapeMask( imFeat, mask, 0.02);

            while outside 

                % Shorten
                for jc = 1 : obj.dim
                    coords = linspace( obj.startPosition(jc), obj.endPosition(jc), 100);
                    obj.endPosition(jc) = coords(95);
                end

                % Check again
                imFeat = obj.simulateFeature( size(mask) );
                outside = Methods.CheckEscapeMask( imFeat, mask);
                obj.GetLength();

            end
            
            
        end
        % }}}
        
        % GetStructInfo {{{
        function feat = GetStructInfo (obj)
            feat.type = obj.type;
            feat.startPosition = obj.startPosition;
            feat.endPosition = obj.endPosition;
            feat.amplitude = obj.amplitude;
            feat.sigma = obj.sigma;
            feat.length = obj.GetLength();
            feat.orientation = obj.GetOrientation();
            feat.ID = obj.ID;
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
            X( X < 1) = 1; X( X > sizeImage(2) ) = sizeImage(2);
            Y( Y < 1) = 1; Y( Y > sizeImage(1) ) = sizeImage(1);
            if dim == 3
                Z = round( linspace( startPoint(3), endPoint(3), round(len) ) );
                Z( Z < 1) = 1; Z( Z > sizeImage(3) ) = sizeImage(3);
            end
            
            if dim == 2
                idxVoxels = sub2ind( sizeImage, Y, X);
            elseif dim ==3
                idxVoxels = sub2ind( sizeImage, Y, X, Z);
            end
            imLine( idxVoxels ) = 1;
            
            % dilate the line with a sphere of large size
            % Assume z-direction is limited (~10 pixels)
            imLine = imdilate( max( imLine, [], 3), strel('disk', radDilate) );
            if dim == 3
                imLine = repmat( imLine, 1, 1, sizeImage(3) );
            end
%             imLine = imdilate( imLine, strel( 'sphere', radDilate) );
            
            % return voxel indices for the dilated line
            lineVox.idx = find( imLine(:) );
            if dim == 2
                [lineVox.y, lineVox.x] = ind2sub( sizeImage, lineVox.idx);
            elseif dim == 3
                [lineVox.y, lineVox.x, lineVox.z] = ind2sub( sizeImage, lineVox.idx);
            end
            
            
        end
        % }}}

    end
end
