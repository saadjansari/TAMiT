classdef Curve < BasicElement

    properties
        startPosition
        cX % coefficients X
        cY % coefficients Y
        cZ % coefficients Z
        order % polynomial order
        length % filament length
        orientation % average orientation
    end

    methods
       
        % Curve {{{
        function obj = Curve( startPosition, coeffs, amplitude, sigma, dim, props2Fit, display)
        % Line : this is the constructor function for a Line. This could be a Microtubule

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( startPosition) || dim ~= length( coeffs) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end
            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'Curve');

            obj.startPosition = startPosition;
            obj.cX = coeffs{1}(1:end-1);
            obj.cY = coeffs{2}(1:end-1);
            obj.order = [ length( obj.cX) , length( obj.cY)];
            if dim == 3
                obj.cZ = coeffs{3}(1:end-1);
                obj.order = [ length( obj.cX), length( obj.cY), length( obj.cZ)];
            end

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props2get)

            % sample props2get
            if nargin==1
                props2get = obj.props2Fit;
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

            if obj.dim == 2
                props2find = {'startPosition', 'cX','cY', 'amplitude', 'sigma'};
            elseif obj.dim == 3
                props2find = {'startPosition', 'cX','cY','cZ', 'amplitude', 'sigma'};
            end
            
            % find the index of start positions
            for jProp = 1 : length( props2find)
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );
                
                % Checking
                if isempty( idxProp)
                    continue;
                end
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVec: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{ jProp} ) = vec( idxProp);

            end
            
%             if any(obj.sigma > 2.01)
%                 disp('curve ub sigma issues')
%                 stoph = 1;
%             end

        end
        % }}}
        
        % simulateFeature {{{
        function [imageOut,error_code] = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: input needed for size of image to simulate the feature in ')
            end
            imageOut = zeros( sizeImage);

            coeffs = obj.GetCoeffFull();
            
            % Simulate a gaussian curve 
            if isfield( obj.params, 'idx')
                if obj.dim == 2
                    imGraph = Cell.drawGaussianCurve2D( coeffs, obj.sigma, imageOut, obj.params.idx, obj.params.x, obj.params.y);
                    error_code = 0;
                elseif obj.dim == 3
                    [imGraph, error_code] = Cell.drawGaussianCurve3D( coeffs, obj.sigma, imageOut, obj.params.idx, obj.params.x, obj.params.y, obj.params.z);
                end
            else
                if obj.dim == 2
                    imGraph = Cell.drawGaussianCurve2D( coeffs, obj.sigma, imageOut);
                    error_code = 0;
                elseif obj.dim == 3
                    [imGraph, error_code] = Cell.drawGaussianCurve3D( coeffs, obj.sigma, imageOut);
                end
            end
            imageFeat = obj.amplitude * mat2gray( imGraph);
            obj.imageSim = imageFeat;
            imageOut = imageFeat + imageOut;

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                error('displayFeature: must provide axes handle to display the feature in')
            end

            % Get (x,y) coordinates of the curve
            coords = obj.GetCoords();

            % Create the curve to display
            line( coords(1,:), coords(2,:), obj.display{:} )

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

            % Get (x,y,z) coordinates of the curve
            coords = obj.GetCoords();

            % Create the curve to display
            line( coords(3,:), coords(1,:), obj.display{:} )

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

            % Get (x,y,z) coordinates of the curve
            coords = obj.GetCoords();

            % Create the curve to display
            line( coords(3,:), coords(2,:), obj.display{:} )

        end
        % }}}
        
        % fillParams {{{
        function obj = fillParams( obj, sizeImage)
            
            coeffs = obj.GetCoeffFull();
            obj.params = Curve.findVoxelsNearCurve( coeffs, sizeImage, 7);
            
        end
        % }}}
        
        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.props2Fit = obj.props2Fit;
            S.startPosition = obj.startPosition;
            S.cX = obj.cX;
            S.cY = obj.cY;
            if obj.dim == 3
                S.cZ = obj.cZ;
            end
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.display = obj.display;

        end
        % }}}
        
        % Get {{{
        
        % GetLength {{{
        function len = GetLength( obj)
            
            coords = obj.GetCoords();
            if obj.dim == 2
                len = sum( sqrt( diff( coords(1,:)).^2 + diff( coords(2,:)).^2 ) );
            elseif obj.dim == 3
                len = sum( sqrt( diff( coords(1,:)).^2 + diff( coords(2,:)).^2 + diff( coords(3,:)).^2 ) );
            end
            obj.length = len;
            
        end
        % }}}
        % GetOrientation {{{
        function orientationXY = GetOrientation( obj)

            coords = obj.GetCoords(); 
            orientationXY = mean( atan( diff( coords(2,:) ) ./ diff( coords(1,:) ) ) );
            %orientationLineZ = obj.fitCoef.z(end-1);

        end
        % }}}
        % GetCoords {{{
        function coords = GetCoords( obj, tv)
            
            coeffs = obj.GetCoeffFull();

            if nargin < 2
                tv = linspace(0,1);
            end

            x = polyval( coeffs{1}, tv )';
            y = polyval( coeffs{2}, tv )';
            if obj.dim == 2
                coords = [x' ; y'];
            elseif obj.dim == 3
                z = polyval( coeffs{3}, tv )';
                coords = [x' ; y'; z'];
            end
            
        end
        % }}}
        % GetCoeffFull {{{
        function coeffs = GetCoeffFull( obj)
            
            if obj.dim == 2
                coeffs = { [obj.cX, obj.startPosition(1)], [ obj.cY, obj.startPosition(2)] };
            elseif obj.dim == 3
                coeffs = { [obj.cX, obj.startPosition(1)], [ obj.cY, obj.startPosition(2)] , [obj.cZ, obj.startPosition(3)]};
            end
            
        end
        % }}}
        % }}}
        
        % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
            % Get 2D projection of feature
            
            % Check object dimensionality
            if obj.dim == 2
                warning('object dimensionality is already 2')
            end
            
            obj.startPosition = obj.startPosition(1:2);
            obj.cZ = [];
            obj.props2Fit = {'cX', 'cY', 'amplitude', 'sigma'};
            
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
            obj.cZ = 0;
            obj.props2Fit = {'cX', 'cY', 'cZ', 'amplitude', 'sigma'};
            
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            
            obj.startPosition(1:2) = obj2D.startPosition(1:2);
            obj.sigma(1:2) = obj2D.sigma(1:2);
            obj.amplitude = obj2D.amplitude;
            obj.cX = obj2D.cX;
            obj.cY = obj2D.cY;
            
        end
        % }}}
        
        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
            % Force features to lie inside the mask. This will
            % shorten any curves whose mtoc is within the mask.

            % Simulate image
            [imFeat,outsideZ] = obj.simulateFeature( size(mask) );
            
            % if escapes from mask, shorten until inside mask
            outsideXY = Methods.CheckEscapeMask( imFeat, mask, 0.02);

            while outsideXY || outsideZ

                % Shorten
                coords = obj.GetCoords( linspace(0,0.95)); 
                coeffs = Curve.estimatePolyCoefficients( coords, obj.order);
                obj.cX = coeffs{1}(1:end-1);
                obj.cY = coeffs{2}(1:end-1);
                if obj.dim == 3
                    obj.cZ = coeffs{3}(1:end-1);
                end

                % Check again
                [imFeat,outsideZ] = obj.simulateFeature( size(mask) );
                outsideXY = Methods.CheckEscapeMask( imFeat, mask);
                obj.GetLength();

            end
            
            % Also ensure that feature is within the z-planes if 3D
            if obj.dim == 3
                if obj.startPosition(3) >= 7
                    obj.cZ(end) = -0.1;
                end
                if obj.startPosition(3) <= 1
                    obj.cZ(end) = 0.1;
                end
            end
            
        end
        % }}}

        % GetStructInfo {{{
        function feat = GetStructInfo(obj)
            feat.type = type;
            feat.startPosition = obj.startPosition;
            feat.amplitude = obj.amplitude;
            feat.sigma = obj.sigma;
            feat.length = obj.GetLength();
            feat.orientation = obj.GetOrientation();
            coords = obj.GetCoords(); 
            feat.endPosition = coords(:,end);
        end
        % }}}

    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Curve')
                error('incorrect type')
            end
            
            if S.dim == 2
                cf = {S.cX, S.cY};
                cf{1}(1+end) = S.startPosition(1);
                cf{2}(1+end) = S.startPosition(2);
            elseif S.dim==3
                cf = {S.cX,S.cY,S.cZ};
                cf{1}(1+end) = S.startPosition(1);
                cf{2}(1+end) = S.startPosition(2);
                cf{3}(1+end) = S.startPosition(3);
            end
            
            obj = Curve( S.startPosition, cf, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);

        end
        % }}}
        
        % findVoxelsNearCurve {{{
        function lineVox = findVoxelsNearCurve( coeffs, sizeImage, rad)
            
            dim = length( sizeImage);
            if length(coeffs) ~= dim
                error('Coeff dim doesnt match dim of image')
            end
            % Create coordinates for numerical integration
            % create parametric coord
            Tvec = linspace(0,1, round(100/5) );Tvec( end) = [];

            % Now for each value of the parameter, we'll do a gauss
            % quadrature numerical integration 
            DeltaT = median( diff(Tvec) );

            % get the offsets for gauss quadrature
            poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
            poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;
            Tsort = sort([  Tvec + poff1,  Tvec + poff2] );

            % Coefficients
            xLoc = polyval( coeffs{1}, Tsort )';
            yLoc = polyval( coeffs{2}, Tsort )';

            % re-parametrize by length
            dParam = 1;
            for jk = 1 : length( xLoc) -1
                dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2)]; 
            end
            dParamNew = linspace( dParam(1), dParam(end), length(dParam) );

            % fit once again, and find now parametric coordinates
            fitX = fit( dParam', xLoc, 'linearinterp'); xDat = fitX( dParamNew);
            fitY = fit( dParam', yLoc, 'linearinterp'); yDat = fitY( dParamNew);
            
            rads = -rad:rad;
            xp = []; yp = []; xp2 = []; yp2 = [];
            for jr = rads
                xp = [ xp ; round(xDat+jr) ];
                yp = [ yp; round(yDat)];
            end
            for jr = rads
                yp2 = [ yp2 ; round(yp+jr) ];
                xp2 = [ xp2 ; round(xp) ];
            end
            pts = unique( [xp2 , yp2] ,'rows');
            lineVox.x = pts(:,1); lineVox.y = pts(:,2);
            
            if dim == 2
                lineVox.idx = sub2ind( sizeImage, pts(:,2), pts(:,1) );
            elseif dim == 3
                ps = zeros( sizeImage(3)*length(pts), 3);
                for zz = 1 : sizeImage(3)
                    ps( 1+length(pts)*(zz-1): zz*length(pts),:) = [pts , zz*ones(length(pts), 1)];
                end
                pts = ps;
                lineVox.z = pts(:,3);
                lineVox.x = pts(:,1); lineVox.y = pts(:,2);
                lineVox.idx = sub2ind( sizeImage, pts(:,2), pts(:,1), pts(:,3) );
            end
            
        end
        % }}}

        % estimatePolyCoefficients {{{
        function coeffs = estimatePolyCoefficients( coords, order)
            % coords: D x N array (D: dimensionality, N: number of points
            % along curve)
            % order: This is the order of the polynomial: positive integer
            
            if length( order) == 1
                order = order * ones( 1, size(coords,1) );
            end
            
            t = linspace(0,1);
            % interpolate for better accuracy
            xi = interp1( linspace(0,1,size(coords,2) ), coords(1,:), t);
            yi = interp1( linspace(0,1,size(coords,2) ), coords(2,:), t);
            
            % fit polynomial of order polyOrder
            if size(coords,1) == 2
                coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) };
            elseif size(coords,1) == 3
                zi = interp1( linspace(0,1, size(coords,2) ), coords(3,:), t);
                coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) , polyfit( t, zi, order(3) ) };
            end

        end
        % }}}
        
        % findCurve {{{
        function coords = findCurve( img, startPoint, orientation)
       
            % Define some search parameters
            stepSize = 5;
            fieldOfView = 40;
            visibility = 10;

            % Use a method to find Curve coordinates given a start point and an orientation.
            coords = Methods.estimateCurveCoords( startPoint(1:2)', orientation, max(img,[],3), stepSize, visibility, fieldOfView,1);

            disp('stop here')
        
        end
        % }}}

    end
end
