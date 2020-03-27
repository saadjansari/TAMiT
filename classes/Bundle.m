classdef Bundle < BasicElement 

    properties
        cX % coefficients X
        cY % coefficients Y
        cZ % coefficients Z
        order % polynomial order
        length % filament length
        orientation % average orientation
        T = [0.3 0.7] % overlap start point, overlap end point. Curve is normalized in length between 0 and 1
        ef = 2 % enhancement factor for amplitude of overlap region
        overlapType = 'Curve'
    end

    methods

        % Bundle {{{
        function obj = Bundle( coeffs, amplitude, sigma, dim, props2Fit, display, t, ef)
        % Bundle : this is the constructor function for an interphase bundle. 

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( coeffs) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end

            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'Bundle');

            obj.cX = coeffs{1};
            obj.cY = coeffs{2};
            obj.order = [ length( obj.cX)-1 , length( obj.cY)-1];
            if dim == 3
                obj.cZ = coeffs{3};
                obj.order = [ length( obj.cX)-1, length( obj.cY)-1, length( obj.cZ)-1];
            end
            obj.T = t;
            obj.ef = ef;

            if length(obj.T) == 1
                obj.overlapType = 'Spot';
            elseif length(obj.T) == 2
                obj.overlapType = 'Curve';
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
                props2find = {'cX','cY', 'amplitude', 'sigma', 'T', 'ef'};
            elseif obj.dim == 3
                props2find = {'cX','cY','cZ', 'amplitude', 'sigma', 'T', 'ef'};
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
                if any(obj.( props2find{ jProp} ) ~= vec( idxProp))
                    stoph = 1;
                end
                obj.( props2find{ jProp} ) = vec( idxProp);
                

            end

        end
        % }}}
        
        % simulateFeature {{{
        function [imageOut,error_code,err_amt] = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: input needed for size of image to simulate the feature in ')
            end
            imageOut = zeros( sizeImage);

            % Simulate a gaussian curve 
            if ~isfield( obj.params, 'idx')
                obj.fillParams(sizeImage);
            end

            % Simulate the three curves
            imageOut = zeros( sizeImage);
            
            % Check if t parameter is consistent, otherwise return with
            % error
            %if any( diff(obj.T) <= 0.05)
                %error_code = 1;
                %return
            %end
            %
            T = sort( obj.T);

            % Get common parameters
            if obj.dim == 2
                spotVars = {obj.params.idx, obj.params.x, obj.params.y};
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y};
                cf = {obj.cX, obj.cY};
            elseif obj.dim == 3
                spotVars = {obj.params.idx, obj.params.x, obj.params.y, obj.params.z};
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y, 'Z', obj.params.z};
                cf = {obj.cX, obj.cY, obj.cZ};
            end

            % Simulate 1st region of curve
            ftype = ['Curve', num2str(obj.dim)];
            [imGraph1, err_code1, err1] = DrawGaussian( obj.sigma, imageOut, ftype, 'Coeff', cf, 'T', [0, T(1)], graphVars{:});
            imageOut = imageOut + obj.amplitude*mat2gray(imGraph1);
            
            % Simulate overlap region 
            ftype = [obj.overlapType, num2str(obj.dim)];
            if strcmp( obj.overlapType, 'Spot')
                try
                    pt = [ polyval( obj.cX, T(1)), polyval( obj.cY, T(1)), polyval( obj.cZ, T(1))];
                catch
                    pt = [ polyval( obj.cX, T(1)), polyval( obj.cY, T(1))];
                end
                imGraph2 = DrawGaussian( obj.sigma, imageIn, ftype, 'Pos', pt, spotVars{:});
            elseif strcmp( obj.overlapType, 'Curve')
                [imGraph2, err_code2, err2] = DrawGaussian( obj.sigma, 0*imageOut, ftype, 'Coeff', cf, 'T', [T(1), T(2)], graphVars{:});
            end
            imageOut = imageOut + obj.ef*obj.amplitude*mat2gray(imGraph2);
            
            % Simulate 3rd region of curve
            [imGraph3, err_code3, err3] = DrawGaussian( obj.sigma, 0*imageOut, ftype, 'Coeff', cf, 'T', [T(2), 1], graphVars{:});
            imageOut = imageOut + obj.amplitude*mat2gray(imGraph3);

            error_code = max( [err_code1, err_code2, err_code3]);
            err_amt = mean( [err1,err2,err3]);
            if err_amt > 0
                stoph = 1;
            end

            obj.imageSim = imageOut;
            if error_code
                stoph=1;
                disp('dang')
            end

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

        % fillParams {{{
        function obj = fillParams( obj, sizeImage)
            if obj.dim == 2
                coeffs = {obj.cX, obj.cY};
            elseif obj.dim == 3
                coeffs = {obj.cX, obj.cY, obj.cZ};
            end
            obj.params = Bundle.findVoxelsNearCurve( coeffs, [0 1], sizeImage, 7);
        end
        % }}}

        % resetT {{{
        function obj = resetT( obj)
           
            % Reset the coefficients to allow the t parameters to lie between 0 and 1.
            t_new = linspace(0,1);
            
            % Get coordinates of whole bundle
            coords = obj.GetCoords();
            
            % Get new coefficients
            coeffs = Bundle.estimatePolyCoefficients( coords, obj.order, linspace(0,1) );
            
            if obj.dim == 2
                error('2d not set up')
            end
            
            if obj.overlapType == 'Curve'
                % Get start/end coordinates of overlap region
                c1 = [ polyval( obj.cX, obj.T(2) ), polyval( obj.cY, obj.T(2) ), polyval( obj.cZ, obj.T(2) )];
                c2 = [ polyval( obj.cX, obj.T(3) ), polyval( obj.cY, obj.T(3) ), polyval( obj.cZ, obj.T(3) )];
                
                coords = [polyval( obj.cX, t_new ) ; polyval( obj.cY, linspace(0,1) ); polyval( obj.cZ, linspace(0,1) )];
                [~, t1] = min( sum( ( repmat( c1, 100, 1)' - coords).^2, 1) ); t1 = 0.01*t1;
                [~, t2] = min( sum( ( repmat( c2, 100, 1)' - coords).^2, 1) ); t2 = 0.01*t2;
                
                obj.T = [t_new(1), t1, t2, t_new(end)];
                
            elseif obj.overlapType == 'Spot'
                % Get coordinates of overlap spot
                c1 = [ polyval( obj.cX, obj.T(2) ), polyval( obj.cY, obj.T(2) ), polyval( obj.cZ, obj.T(2) )];
                
                coords = [polyval( obj.cX, t_new ) ; polyval( obj.cY, linspace(0,1) ); polyval( obj.cZ, linspace(0,1) )];
                [~, t1] = min( sum( ( repmat( c1, 100, 1)' - coords).^2, 1) ); t1 = 0.01*t1;
                
                obj.T = [t_new(1), t1, t_new(end)];
            
            end
                
            obj.cX = coeffs{1};
            obj.cY = coeffs{2};
            if obj.dim == 3
                obj.cZ = coeffs{3};
            end
            
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
            
            if nargin < 2
                tv = linspace( 0, 1);
            end
            x = polyval( obj.cX, tv )';
            y = polyval( obj.cY, tv )';
            if obj.dim == 2
                coords = [x' ; y'];
            elseif obj.dim == 3
                z = polyval( obj.cZ, tv )';
                coords = [x' ; y'; z'];
            end
            
        end
        % }}}
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
                coords = obj.GetCoords( linspace(0.02,0.98)); 
                coeffs = Bundle.estimatePolyCoefficients( coords, obj.order, linspace(0, 1) );
                obj.cX = coeffs{1};
                obj.cY = coeffs{2};
                if obj.dim == 3
                    obj.cZ = coeffs{3};
                end

                % Check again
                [imFeat,outsideZ] = obj.simulateFeature( size(mask) );
                outsideXY = Methods.CheckEscapeMask( imFeat, mask);
                obj.GetLength();

            end
            
            % Also ensure that feature is within the z-planes if 3D
            if obj.dim == 3
                if obj.cZ(end) >= 7
                    obj.cZ(end) = 7;
                    obj.cZ(end-1) = -0.3;
                end
                if obj.cZ(end) <= 1
                    obj.cZ(end) = 1;
                    obj.cZ(end-1) = +0.3;
                end
            end
            
        end
        % }}}

        % GetStructInfo {{{
        function feat = GetStructInfo(obj)
            feat.type = type;
            feat.cX = obj.cX;
            feat.cY = obj.cY;
            feat.cZ = obj.cZ;
            feat.amplitude = obj.amplitude;
            feat.sigma = obj.sigma;
            feat.length = obj.GetLength();
            feat.orientation = obj.GetOrientation();
            coords = obj.GetCoords(); 
            feat.endPosition = coords(:,end);
        end
        % }}}

        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.props2Fit = obj.props2Fit;
            S.cX = obj.cX;
            S.cY = obj.cY;
            if obj.dim == 3
                S.cZ = obj.cZ;
            end
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.t = obj.T;
            S.ef = obj.ef;
            S.display = obj.display;
            S.overlapType = obj.overlapType;

        end
        % }}}

    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Bundle')
                error('incorrect type')
            end
            
            if S.dim == 2
                cf = {S.cX, S.cY};
            elseif S.dim==3
                cf = {S.cX,S.cY,S.cZ};
            end
            
            obj = Bundle( cf, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display, S.t, S.ef);

        end
        % }}}
        
        % findVoxelsNearCurve {{{
        function lineVox = findVoxelsNearCurve( coeffs, t, sizeImage, rad)
            
            dim = length( sizeImage);
            if length(coeffs) ~= dim
                error('Coeff dim doesnt match dim of image')
            end
            % Create coordinates for numerical integration
            % create parametric coord
            Tvec = linspace( t(1), t(2), round(100/5) );Tvec( end) = [];

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
                try
                lineVox.idx = sub2ind( sizeImage, pts(:,2), pts(:,1), pts(:,3) );
                catch
                    stoph = 1;
                end
            end
            
        end
        % }}}

        % estimatePolyCoefficients {{{
        function coeffs = estimatePolyCoefficients( coords, order, t)
            % coords: D x N array (D: dimensionality, N: number of points
            % along curve)
            % order: This is the order of the polynomial: positive integer
            
            if length( order) == 1
                order = order * ones( 1, size(coords,1) );
            end
            if nargin < 3
                t = linspace(0,1);
            end

            % interpolate for better accuracy
            xi = interp1( linspace(0,1,size(coords,2) ), coords(1,:), t, 'linear', 'extrap');
            yi = interp1( linspace(0,1,size(coords,2) ), coords(2,:), t, 'linear', 'extrap');
            
            % fit polynomial of order polyOrder
            if size(coords,1) == 2
                coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) };
            elseif size(coords,1) == 3
                zi = interp1( linspace(0,1, size(coords,2) ), coords(3,:), t, 'linear', 'extrap');
                coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) , polyfit( t, zi, order(3) ) };
            end

        end
        % }}}

    end

end

