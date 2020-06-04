classdef BundleNew < BasicElement 
    % This is a bundle which represents a curve that is linear in the Z
    % dimension. The XY dimensions are modeled with:
    % 1) Start point - 3 params
    % 2) Initial unit tangent vector - 3 params
    % 3) Normal vector as a function of time (first order polynomial) - 2x2
    % = 4 params (for both sides)
    % 4) Maximum time - 1x2 = 2 params (gets at length of microtubule)
    % 5) Overlap end time (from center of overlap zone) - 1 param
    % 6) Enhancement factor - 1 param
    % 7) Sigma - 3 params (why not fix this to be PSF?)
    % 8) amplitude - 1 param
    % Total parameters = 3+3+4+2+1+1+3 = 18 params

    properties
        origin % origin of bundle (center of overlap zone)
        thetaInit % initial tangent theta and phi
        normalVec % 2x2 array (two rows represent the normal vector polynomial coeffs for the two extensions of the bundle.
        L = [8 8] % length of both extensions.
        T = 1 % overlap end point from origin. Curve is normalized in length between 0 and 1
        ef = 2 % enhancement factor for amplitude of overlap region
        length % filament length
        orientation % average orientation
        two_sided = 1
    end

    methods

        % Bundle {{{
        function obj = BundleNew( origin, thetaInit, normalVec, T, L, ef, amplitude, sigma, dim, props2Fit, display)
        % Bundle : this is the constructor function for an interphase bundle. 

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( origin) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input origin')
            end

            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'BundleNew');

            obj.origin = origin;
%             obj.tanInit = [ tanInit(1:2) ; tanInit( 3:4)];
%             obj.normalVec = [ normalVec(1:end/2) ; normalVec( end/2:end)];
            obj.thetaInit = thetaInit;
            obj.normalVec = normalVec;
            obj.L = L; obj.L( obj.L < 8) = 9;
            obj.T = T;
            obj.ef = ef;
            if length(obj.L) == 1
                obj.two_sided = 0; % only 1 extension
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

            props2find = {'origin','tanInit', 'normalVec', 'T', 'L', 'ef', 'amplitude', 'sigma'};
            
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

            % Simulate a gaussian curve 
            if ~isfield( obj.params, 'idx')
                obj.fillParams(sizeImage);
            end

            % Simulate the three curves
            imageOut = zeros( sizeImage);

            % Get common parameters
            if obj.dim == 2
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y};
            elseif obj.dim == 3
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y, 'Z', obj.params.z};
            end
            errs = []; errcodes = [];
            ftype = ['Curve', num2str(obj.dim) 'Coords'];
            % Simulate 1st region of curve
            if obj.L(1)>obj.T
                [imGraph1, ec1, err1] = DrawGaussian( obj.sigma, imageOut, ftype, 'Coord', obj.GetCoords(1,obj.L(1),obj.T), graphVars{:});
                imageOut = imageOut + obj.amplitude*mat2gray(imGraph1);
                errs = [errs, err1]; errcodes = [errcodes, ec1];
            end
            % Simulate overlap region 1
            [imGraph2, ec2, err2] = DrawGaussian( obj.sigma, 0*imageOut, ftype, 'Coord', obj.GetCoords(1,obj.T,0), graphVars{:});
            imageOut = imageOut + obj.ef*obj.amplitude*mat2gray(imGraph2);
            errs = [errs, err2]; errcodes = [errcodes, ec2];
            
            if obj.two_sided
                % Simulate overlap region 2 
                [imGraph3, ec3, err3] = DrawGaussian( obj.sigma, 0*imageOut, ftype, 'Coord', obj.GetCoords(2,obj.T,0), graphVars{:});
                imageOut = imageOut + obj.ef*obj.amplitude*mat2gray(imGraph3);
                errs = [errs, err3]; errcodes = [errcodes, ec3];

                % Simulate 3rd region of curve
                if obj.L(2) > obj.T
                    [imGraph4, ec4, err4] = DrawGaussian( obj.sigma, 0*imageOut, ftype, 'Coord', obj.GetCoords(2,obj.L(2),obj.T), graphVars{:});
                    imageOut = imageOut + obj.amplitude*mat2gray(imGraph4);
                    errs = [errs, err4]; errcodes = [errcodes, ec4];
                end
            end
            error_code = max( errcodes);
            err_amt = mean( errs);
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

            cx = []; cy = [];
            if obj.L(1)>obj.T
                c1 = obj.GetCoords(1,obj.L(1),obj.T);
                cx = [cx, c1(1,end:-1:1)]; cy = [cy, c1(2,end:-1:1)];
            end
            c2 = obj.GetCoords(1,obj.T,0);
            cx = [cx, c2(1,end:-1:1)]; cy = [cy, c2(2,end:-1:1)];

            if obj.two_sided
                c3 = obj.GetCoords(2,obj.T,0);
                cx = [cx, c3(1,:)]; cy = [cy, c3(2,:)];
                if obj.L(2) > obj.T
                    c4 = obj.GetCoords(2,obj.L(2),obj.T);
                    cx = [cx, c4(1,:)]; cy = [cy, c4(2,:)];
                end
            end
            coords = [cx; cy];
            % Get (x,y) coordinates of the curve
            %coords = obj.GetCoords();

            % Create the curve to display
            line( coords(1,:), coords(2,:), obj.display{:} )

        end
        % }}}

        % fillParams {{{
        function obj = fillParams( obj, sizeImage)

            obj.params = BundleNew.findVoxelsNearCurve( obj.GetCoords(), sizeImage, 10);
            
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
        function cc = GetCoords( obj, whichExt, tmax, tmin)
            % whichExt = 0 (both extensions), 1 (1st extension), 2 (2nd
            % extension)
            if nargin < 2
                whichExt = 0;
            end
            if obj.two_sided == 0
                whichExt = 1;
            end
            if nargin < 3
                tmax = obj.L;
            end
            
            if nargin < 4
                tmin = 0;
            end
            
            Rot = [ 0 -1; 1 0]; % Rotation matrix
            switch whichExt
                case 0
                    nV = [obj.normalVec(1:2); obj.normalVec(3:4)];
                    thInit = [ obj.thetaInit(1:2); obj.thetaInit(3:4)];
                    if nargin < 3
                        tmax = obj.L;
                    end
                    if length(tmax) == 1
                        tmax(2) = tmax(1);
                    end
                    % For each extension
                    xxs = {};
                    for jj = 1:2
                        
                        tanVec = [cos(thInit(jj,1)), sin(thInit(jj,1)), cos( thInit(jj,2))];
                        % construct time vector
                        t = 0:0.1:tmax(jj); [~,idxStart] = min( abs(t-tmin));
                        
                        % Acceleration function discretized in time
                        acc = nV(jj,1) + nV(jj,2)*t;

                        % initialization
                        xx = zeros( obj.dim, length(t) ); xx(:,1) = obj.origin;
                        vv = zeros( 2, length(t) ); vv(:,1) = tanVec(1:2);

                        % iterate
                        for jt = 2 : length(t)
                            vv(:,jt) = vv(:,jt-1) + acc( jt)*0.1*Rot* vv(:,jt-1) / norm( vv(:,jt-1));
                            v_mean = 1/2 * ( vv(:,jt) + vv(:,jt-1) );
                            xx(1:2,jt) = xx(1:2,jt-1) + 0.1*v_mean/norm(v_mean);
                            if obj.dim == 3
                                xx(3,jt) = xx(3,jt-1) + tanVec(3)*0.1;
                            end
                        end
                        xxs{jj} = xx;
                    end
                    cc = [flip( xxs{1}(:,2:end), 2) xxs{2} ];
                    
                case 1
                    nV = obj.normalVec(1:2);
                    thInit = obj.thetaInit(1:2);
                    if nargin < 3
                        tmax = obj.L(1);
                    end
                    
                    tanVec = [cos(thInit(1,1)), sin(thInit(1,1)), cos( thInit(1,2))];
                    % construct time vector
                    t = 0:0.1:tmax(1);  [~,idxStart] = min( abs(t-tmin));

                    % Acceleration function discretized in time
                    acc = nV(1,1) + nV(1,2)*t;

                    % initialization
                    xx = zeros( obj.dim, length(t) ); xx(:,1) = obj.origin;
                    vv = zeros( 2, length(t) ); vv(:,1) = tanVec(1:2);

                    % iterate
                    for jt = 2 : length(t)
                        vv(:,jt) = vv(:,jt-1) + acc( jt)*0.1*Rot* vv(:,jt-1) / norm( vv(:,jt-1));
                        v_mean = 1/2 * ( vv(:,jt) + vv(:,jt-1) );
                        xx(1:2,jt) = xx(1:2,jt-1) + 0.1*v_mean/norm(v_mean);
                        if obj.dim == 3
                            xx(3,jt) = xx(3,jt-1) + tanVec(3)*0.1;
                        end
                    end
                    cc = xx(:,idxStart:end);
                    
                case 2
                    nV = [obj.normalVec(1:2); obj.normalVec(3:4)];
                    thInit = [ obj.thetaInit(1:2); obj.thetaInit(3:4)];
                    if nargin < 3
                        tmax = obj.L(2);
                    end
                    
                    tanVec = [cos(thInit(2,1)), sin(thInit(2,1)), cos( thInit(2,2))];
                    % construct time vector
                    t = tmin:0.1:tmax(1);  [~,idxStart] = min( abs(t-tmin));

                    % Acceleration function discretized in time
                    acc = nV(2,1) + nV(2,2)*t;

                    % initialization
                    xx = zeros( obj.dim, length(t) ); xx(:,1) = obj.origin;
                    vv = zeros( 2, length(t) ); vv(:,1) = tanVec(1:2);

                    % iterate
                    for jt = 2 : length(t)
                        vv(:,jt) = vv(:,jt-1) + acc( jt)*0.1*Rot* vv(:,jt-1) / norm( vv(:,jt-1));
                        v_mean = 1/2 * ( vv(:,jt) + vv(:,jt-1) );
                        xx(1:2,jt) = xx(1:2,jt-1) + 0.1*v_mean/norm(v_mean);
                        if obj.dim == 3
                            xx(3,jt) = xx(3,jt-1) + tanVec(3)*0.1;
                        end
                    end
                    cc = xx(:,idxStart:end);
                    
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
                obj.L = obj.L-1;
                obj.T = max([2 obj.T-1]);

                % Check again
                [imFeat,outsideZ] = obj.simulateFeature( size(mask) );
                outsideXY = Methods.CheckEscapeMask( imFeat, mask);
                obj.GetLength();

            end
            
%             % Also ensure that feature is within the z-planes if 3D
%             if obj.dim == 3
%                 if obj.cZ(end) >= 7
%                     obj.cZ(end) = 7;
%                     obj.cZ(end-1) = -0.3;
%                 end
%                 if obj.cZ(end) <= 1
%                     obj.cZ(end) = 1;
%                     obj.cZ(end-1) = +0.3;
%                 end
%             end
            
        end
        % }}}

        % GetStructInfo {{{
        function feat = GetStructInfo(obj)
            feat.type = type;
            feat.origin = obj.origin;
            feat.tanInit = obj.tanInit;
            feat.normalVec = obj.normalVec;
            feat.L = obj.L;
            feat.T = obj.T;
            feat.ef = obj.ef;
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
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.origin = obj.origin;
            S.thetaInit = obj.thetaInit;
            S.normalVec = obj.normalVec;
            S.L = obj.L;
            S.T = obj.T;
            S.ef = obj.ef;
            S.display = obj.display;

        end
        % }}}

    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'BundleNew')
                error('incorrect type')
            end          
            
            obj = BundleNew( S.origin, S.thetaInit, S.normalVec, S.T, S.L, S.ef, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);

        end
        % }}}
        
        % findVoxelsNearCurve {{{
        function lineVox = findVoxelsNearCurve( coords, sizeImage, rad)
            
            dim = length( sizeImage);
            if size(coords,1) ~= dim
                error('Coord size doesnt match dim of image')
            end
            
            rads = -rad:rad;
            xp = []; yp = []; xp2 = []; yp2 = [];
            for jr = rads
                xp = [ xp ; round(coords(1,:)'+jr) ];
                yp = [ yp; round(coords(2,:)')];
            end
            for jr = rads
                yp2 = [ yp2 ; round(coords(2,:)'+jr) ];
                xp2 = [ xp2 ; round(coords(1,:)') ];
            end
            % Add half circle near endpoints
            imt = zeros(sizeImage(1:2)); 
            imt(round(coords(2,1)),round(coords(1,1)))=1; imt(round(coords(2,end)),round(coords(1,end)))=1; 
            imt = imdilate( imt, strel('disk',rad)); idxEnd = find(imt);
            [ye, xe] = ind2sub( sizeImage(1:2), idxEnd);
            xp2 = [xe; xp2]; yp2 = [ye; yp2];
            
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

    end

end

