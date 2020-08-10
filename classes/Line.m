classdef Line < BasicElement

    properties
        startPosition
        endPosition
        length
        theta % [phi,theta] physics convention (phi:0-2pi, theta:0-pi)
        repr = 'cartesian' % cartesian or spherical
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
            if obj.length < 4
                warning('line length is less than 4')
            end
            
            % find orientations
            obj.theta = atan2( obj.endPosition(2)-obj.startPosition(2),obj.endPosition(1)-obj.startPosition(1) );
            if obj.dim == 3
                phi = acos( (obj.endPosition(3)-obj.startPosition(3))/obj.length );
                obj.theta = [obj.theta, phi];
            end
            obj.SetBounds();
            
        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels, ub,lb] = getVec( obj, props2get)

            % sample props2get
            if nargin==1
                if strcmp( obj.repr, 'cartesian')
                    props2get = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
                elseif strcmp( obj.repr, 'spherical')
                    props2get = {'startPosition','length', 'theta','amplitude', 'sigma'};
                end
            end
            
            % make sure props input matches properties defined in class
            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
            end

            % Get vector of Properties
            vec = []; ub = []; lb = [];
            for jProp = 1 : length( props2get)
                if strcmp( obj.repr, 'cartesian') && ( strcmp(props2get{jProp},'theta') || strcmp(props2get{jProp},'length'))
                    continue
                elseif strcmp( obj.repr, 'spherical') && strcmp(props2get{jProp},'endPosition')
                    continue
                end
                vec = [ vec, obj.( props2get{jProp} ) ];
                try
                    ub = [ ub, obj.bounds.ub.( props2get{jProp} ) ];
                    lb = [ lb, obj.bounds.lb.( props2get{jProp} ) ];
                end
            end

            % Also get a string array with property names
            vecLabels = {};
            for jProp = 1 : length( props2get)
                if numel( obj.( props2get{jProp} ) ) ~= length( obj.( props2get{jProp} ) )
                    error('getVec : property selected is a non-singleton matrix')
                end
                if strcmp( obj.repr, 'cartesian') && ( strcmp(props2get{jProp},'theta') || strcmp(props2get{jProp},'length'))
                    continue
                elseif strcmp( obj.repr, 'spherical') && strcmp(props2get{jProp},'endPosition')
                    continue
                end
                numRep = length( obj.( props2get{jProp} ) );
                labelRep = cell( 1, numRep);
                labelRep(:) = props2get(jProp);
                vecLabels = { vecLabels{:}, labelRep{:} };
            end
            if isempty(ub)
                clearvars ub lb
            end

        end
        % }}}
        
        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)

            props2find = {'startPosition', 'endPosition', 'amplitude', 'sigma','theta','phi', 'length'};

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
            
            % Update props
            if strcmp( obj.repr,'spherical')
                if obj.dim == 3
                    obj.endPosition = obj.startPosition + obj.length* [sin(obj.theta(2))*cos(obj.theta(1)), sin(obj.theta(2))*sin(obj.theta(1)), cos(obj.theta(2))];
                elseif obj.dim == 2
                    obj.endPosition = obj.startPosition + obj.length* [cos(obj.theta(1)), sin(obj.theta(1))];
                end
            end

        end
        % }}}
        
        % simulateFeature {{{
        function [imageOut,ec,err] = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: input needed for size of image to simulate the feature in ')
            end
            imageOut = zeros( sizeImage);

            % Simulate a gaussian line
            if isfield( obj.params, 'idx')
                if obj.dim == 2
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Line2', 'PosStart', obj.startPosition, 'PosEnd', obj.endPosition,'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y);
                    imageOut = obj.amplitude * mat2gray( imGraph);
                elseif obj.dim == 3
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Line3', 'PosStart', obj.startPosition, 'PosEnd', obj.endPosition,'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y, 'Z', obj.params.z);
                    imageOut = obj.amplitude * mat2gray( imGraph);
                end
            else
                if obj.dim == 2
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Line2', 'PosStart', obj.startPosition, 'PosEnd', obj.endPosition);
                    imageOut = obj.amplitude * mat2gray( imGraph);
                elseif obj.dim == 3
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Line3', 'PosStart', obj.startPosition, 'PosEnd', obj.endPosition);
                    imageOut = obj.amplitude * mat2gray( imGraph);
                end
            end
            obj.imageSim = imageOut;
            %imageOut( imageFeat > imageIn) = imageFeat( imageFeat > imageIn);

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, sizeZ)

            if nargin < 2
                error('displayFeature: must provide axes handle to display the feature in')
            end

            % 2D color plot for 3D information
            if obj.dim==3 && nargin==3
                coords = [linspace(obj.startPosition(1), obj.endPosition(1),100); ...
                    linspace(obj.startPosition(2), obj.endPosition(2),100); ...
                    linspace(obj.startPosition(3), obj.endPosition(3),100)];
                cm = cool;
                cols = round((coords(3,:)/sizeZ)*length(cm));
                cols( cols < 1) = 1; cols(cols>256) = 256;
                col = cm( cols, :);
                col = [ permute(col, [3 1 2]); permute(col, [3 1 2])];
                z = zeros([ 1, size( coords,2)]);
                surface([coords(1,:);coords(1,:)],[coords(2,:);coords(2,:)],[z;z],col,...
                        'facecol','no',...
                        'edgecol','interp',...
                        'linew',4);
                colorbar('Ticks',linspace(0,1,sizeZ),'TickLabels',1:sizeZ)
            else
                % Create the line to display
                line( [obj.startPosition(1) obj.endPosition(1)], [obj.startPosition(2) obj.endPosition(2)], obj.display{:} )
            end
            
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
%                 for jc = 1 : obj.dim
%                     coords = linspace( obj.startPosition(jc), obj.endPosition(jc), 100);
%                     obj.endPosition(jc) = coords(95);
%                 end
                obj.absorbVec( [obj.length*0.95], {'length'});

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
            feat.phi = obj.phi;
            feat.theta = obj.theta;
%             feat.orientation = obj.GetOrientation();
            feat.ID = obj.ID;
        end
        % }}}
        
        % SetBounds {{{
        function SetBounds( obj)
           
            % amplitude
            ub.amplitude = 1;
            lb.amplitude = 0;
            
            % sigma % positions % theta
            if obj.dim == 3
                if strcmpi(obj.label, 'spindle')
                	ub.sigma = [3.0 3.0 2.0];
                else
                    ub.sigma = [2.0 2.0 2.0];
                end
                lb.sigma = [1.2 1.2 1.0];
                ub.startPosition = [150 150 7];
                ub.endPosition = [150 150 7];
                lb.startPosition = [1 1 1];
                lb.endPosition = [1 1 1];
                ub.theta = [pi pi];
                lb.theta = [-pi 0];
            elseif obj.dim == 2
                if strcmpi(obj.label, 'spindle')
                	ub.sigma = [3.0 3.0];
                else
                    ub.sigma = [2.0 2.0];
                end
                lb.sigma = [1.2 1.2];
                ub.startPosition = [150 150];
                ub.endPosition = [150 150];
                lb.startPosition = [1 1];
                lb.endPosition = [1 1];
                ub.theta = [pi];
                lb.theta = [-pi];
            end
            
            % length
            ub.length = 300;
            lb.length = 4;
            
            obj.bounds.lb = lb;
            obj.bounds.ub = ub;
            
        end
        % }}}
        
        % preOptimize {{{
        function obj = preOptimize(obj, imOrg, imBkg)
           
            % optimize sigma
            % sx
            res = []; sx = linspace(1.2,2.5,15); sx0 = obj.sigma(1); sy0 = obj.sigma(2);
            for ix = sx
                obj.sigma(1) = ix; obj.sigma(2) = ix;
                imSim = imBkg + obj.simulateFeature( size(imBkg));
                res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
            end
            [~, idx] = min( res); 
            obj.sigma(1) = sx(idx);
            obj.sigma(2) = sx(idx);

            % sz
            if obj.dim == 3
                res = []; sz = linspace(1.0,2.0,10); sz0 = obj.sigma(3);
                for ix = sz
                    obj.sigma(3) = ix;
                    imSim = imBkg + obj.simulateFeature( size(imBkg));
                    res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
                end
                [~, idx] = min( res); obj.sigma(3) = sz(idx);
            end
            %{
             {% optimize amp 
             {res = []; amps = linspace( 0, max( imOrg(:))-imBkg(1), 20); 
             {A0 = obj.amplitude; 
             {for ia = amps
             {    obj.amplitude = ia; 
             {    imSim = imBkg + obj.simulateFeature( size(imBkg));
             {    res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
             {end
             {[~,idA] = min( min(res,[],2) ); obj.amplitude = amps(idA); 
             %}
            
        end
        % }}}

        function SetTheta(obj, theta)
            
            if length(theta) ~= length(obj.theta)
                error('incorrect length of theta')
            end
            
            obj.theta = theta;
            if obj.dim == 3
                obj.endPosition = obj.startPosition + obj.length* [sin(obj.theta(2))*cos(obj.theta(1)), sin(obj.theta(2))*sin(obj.theta(1)), cos(obj.theta(2))];
            elseif obj.dim == 2
                obj.endPosition = obj.startPosition + obj.length* [cos(obj.theta(1)), sin(obj.theta(1))];
            end
                
        end
        function SetLength(obj, len)
            
            if length(len) ~= length(obj.length)
                error('incorrect length of len')
            end
            
            obj.length = len;
            if obj.dim == 3
                obj.endPosition = obj.startPosition + obj.length* [sin(obj.theta(2))*cos(obj.theta(1)), sin(obj.theta(2))*sin(obj.theta(1)), cos(obj.theta(2))];
            elseif obj.dim == 2
                obj.endPosition = obj.startPosition + obj.length* [cos(obj.theta(1)), sin(obj.theta(1))];
            end
                
        end
        function SetEndPosition(obj, endP)
            
            if length(endP) ~= length(obj.endPosition)
                error('incorrect length of endPosition')
            end
            
            obj.endPosition = endP;
            obj.theta = atan2( obj.endPosition(2)-obj.startPosition(2),obj.endPosition(1)-obj.startPosition(1) );
            if obj.dim == 3
                phi = acos( (obj.endPosition(3)-obj.startPosition(3))/obj.length );
                obj.theta = [obj.theta, phi];
            end
            obj.length = norm( obj.startPosition - obj.endPosition);
            
        end
        function SetStartPosition(obj, startP)
            
            if length(startP) ~= length(obj.startPosition)
                error('incorrect length of endPosition')
            end
            
            obj.startPosition = startP;
            obj.theta = atan2( obj.endPosition(2)-obj.startPosition(2),obj.endPosition(1)-obj.startPosition(1) );
            if obj.dim == 3
                phi = acos( (obj.endPosition(3)-obj.startPosition(3))/obj.length );
                obj.theta = [obj.theta, phi];
            end
            obj.length = norm( obj.startPosition - obj.endPosition);
            
        end
        
    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Line')
                error('incorrect type')
            end

            obj = Line( S.startPosition, S.endPosition, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);

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
