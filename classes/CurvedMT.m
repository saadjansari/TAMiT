classdef CurvedMT < BasicElement 
    % This is a curve which represents a curve that is linear in the Z
    % dimension. The XY dimensions are modeled with:
    % 1) Start point - 3 params
    % 2) Initial unit tangent vector - 3 params
    % 3) Normal vector as a function of time (first order polynomial) - 2 params 
    % 4) Maximum time - 1 params (gets at length of microtubule)
    % 7) Sigma - 3 params (why not fix this to be PSF?)
    % 8) amplitude - 1 param
    % Total parameters = 3+3+2+1+3 = 12 params

    properties
        origin % origin of bundle (center of overlap zone)
        thetaInit % initial tangent theta and phi
        curvature % nx1 array
        model
        L = 8 
        length % filament length
        orientation % average orientation
        err_origin % origin of bundle (center of overlap zone)
        err_thetaInit % initial tangent theta and phi
        err_curvature % nx1 array
        err_L
    end

    methods

        % CurvedMT {{{
        function obj = CurvedMT( model, origin, thetaInit, curvature, L, amplitude, sigma, dim, props2Fit, display)
        % Bundle : this is the constructor function for an interphase bundle. 

            % Ensure dim matches image dimensionality and positions dimensionality
            if dim ~= length( origin) || dim ~= length( sigma)
                error( 'Feature: input argument dim does not match dimensionality of input origin')
            end

            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'CurvedMT');
            
            obj.model = model;
            obj.origin = origin;
            obj.thetaInit = thetaInit;
            obj.curvature = curvature;
            obj.L = L; 
            obj.SetBounds();

        end
        % }}}

        % getVec {{{
        function [vec, vecLabels, ub, lb] = getVec( obj, props2get)

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
            if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                temps.origin = obj.origin(3);
                temps.thetaInit = obj.thetaInit(2);
                temps.sigma = obj.sigma(3);
                temps.amplitude = obj.amplitude;
                temps.L = obj.L;
                temps.curvature = obj.curvature;
                temps.bounds.ub.origin = obj.bounds.ub.origin(3);
                temps.bounds.ub.thetaInit = obj.bounds.ub.thetaInit(2);
                temps.bounds.ub.sigma = obj.bounds.ub.sigma(3);
                temps.bounds.ub.amplitude = obj.bounds.ub.amplitude;
                temps.bounds.ub.L = obj.bounds.ub.L;
                temps.bounds.ub.curvature = obj.bounds.ub.curvature;
                temps.bounds.lb.origin = obj.bounds.lb.origin(3);
                temps.bounds.lb.thetaInit = obj.bounds.lb.thetaInit(2);
                temps.bounds.lb.sigma = obj.bounds.lb.sigma(3);
                temps.bounds.lb.amplitude = obj.bounds.lb.amplitude;
                temps.bounds.lb.L = obj.bounds.lb.L;
                temps.bounds.lb.curvature = obj.bounds.lb.curvature;
            end
            vec = []; ub=[]; lb=[];
            for jProp = 1 : length( props2get)
                if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                    vec = [ vec, temps.( props2get{jProp} ) ];
                    ub = [ ub, temps.bounds.ub.( props2get{jProp} ) ];
                    lb = [ lb, temps.bounds.lb.( props2get{jProp} ) ];
                else
                    vec = [ vec, obj.( props2get{jProp} ) ];
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
                if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                    numRep = length( temps.( props2get{jProp} ) );                 
                else
                    numRep = length( obj.( props2get{jProp} ) );
                end
                labelRep = cell( 1, numRep);
                labelRep(:) = props2get(jProp);
                vecLabels = { vecLabels{:}, labelRep{:} };
            end

        end
        % }}}
        
        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels, errBoolean)

            if nargin < 4
                errBoolean = 0;
            end
            
            props2find = {'origin','thetaInit', 'curvature', 'L', 'amplitude', 'sigma'};
            
            % find the index of start positions
            for jProp = 1 : length( props2find)
                
                if errBoolean
                    propCurr = ['err_',props2find{ jProp}];
                else
                    propCurr = props2find{ jProp};
                end
                
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );
                
                % Checking
                if isempty( idxProp)
                    continue;
                end
                if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                    obj.( propCurr ) = obj.( props2find{ jProp} );
                    obj.( propCurr )(1+end-length(idxProp):end) = vec( idxProp);
                    
                else
                    if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                        error( 'absorbVec: length of vector props to absorb does not match the old property size')
                    end

                    % Set final property
                    obj.( propCurr ) = vec( idxProp);
                end

            end
            obj.sigma(1) = mean( obj.sigma(1:2));
            obj.sigma(2) = obj.sigma(1);
            
            % Check if any property is very close but below or above bound,
            % then set it equal to bound
            for jProp = 1 : length( props2find)
                propCurr = props2find{ jProp};
                pc = obj.( propCurr );
                pcu = obj.bounds.ub.(propCurr);
                pcl = obj.bounds.lb.(propCurr);
                for jv = 1: length(pc)
                    if pc(jv) > pcu(jv) && pc(jv) - pcu(jv) < 1e-8
                        pc(jv) = pcu(jv);
                    end
                    if pc(jv) < pcl(jv) && pcl(jv) - pc(jv) < 1e-8
                        pc(jv) = pcl(jv);
                    end
                end
                obj.( propCurr ) = pc;
            end
            
        end
        % }}}
        
        % simulateFeature {{{
        function [imageOut,ec,err] = simulateFeature( obj, sizeImage)
            % Simulate a gaussian curve 

            if nargin < 2
                error('simulateFeature: input needed for size of image to simulate the feature in ')
            end
            if ~isfield( obj.params, 'idx')
                obj.fillParams(sizeImage);
            end

            % Get common parameters
            ftype = ['Curve', num2str(obj.dim) 'Coords'];
            imageOut = zeros( sizeImage);
            if obj.dim == 2
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y};
            elseif obj.dim == 3
                graphVars = {'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y, 'Z', obj.params.z};
            end

            % Simulate curve
            [imSim, ec, err] = DrawGaussian( obj.sigma, imageOut, ftype, 'Coord', obj.GetCoords(), graphVars{:});
            imageOut = obj.amplitude*mat2gray(imSim);
            
            % What to do if there is an error in Z?
            % Options:
            % 1. Increase intensity of the penetrating tip pixel, scaling it with
            %    the error amount.
            
            % Find all z-indices of the max-intensity pixels
            if err > 0
                errorPlane = imageOut(:,:,ec);
                idx = find( errorPlane == max( errorPlane(:) ) );
                [yidx,xidx] = ind2sub( size( errorPlane), idx);
                imageOut( yidx,xidx, ec) = imageOut( yidx,xidx, ec)*(1 +err);
            end
            obj.imageSim = imageOut;

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, sizeZ)

            if nargin < 2
                error('displayFeature: must provide axes handle to display the feature in')
            end

            % Get (x,y) coordinates of the curve
            coords = obj.GetCoords();
            
            
            % 2D color plot for 3D information
            if obj.dim==3 && nargin==3
                
                % Force coords betwwen planes
                coords(3, find(coords(3,:)>sizeZ)) = sizeZ;
                coords(3, find(coords(3,:)<1)) = 1;
                
                cm = hsv;
                col = cm( round((coords(3,:)/sizeZ)*length(cm)), :);
                col = [ permute(col, [3 1 2]); permute(col, [3 1 2])];
                z = zeros([ 1, size( coords,2)]);
                surface([coords(1,:);coords(1,:)],[coords(2,:);coords(2,:)],[z;z],col,...
                        'facecol','no',...
                        'edgecol','interp',...
                        'linew',3);
                % colorbar('Ticks',linspace(0,1,sizeZ),'TickLabels',1:sizeZ)
                
            else
                % Create the curve to display
                line( coords(1,:), coords(2,:), obj.display{:} )
            end

        end
        % }}}
        
        % displayFeature {{{
        function ax = displayFeature_color_by_curvature( obj, ax)

            if nargin < 2
                error('displayFeature: must provide axes handle to display the feature in')
            end

            % Get (x,y) coordinates of the curve
            coords = obj.GetCoords();
            
            % 2D color plot for 3D information
            % Extract local curvature
            t = linspace(0, obj.L, size(coords,2));
            acc = obj.MakeModelFunction();
            K = abs( acc(t));
            % range of K from 0 to 0.05
            K(K>0.15) = 0.15;
                
            cm = brighten(jet,0.4);
            vals = round( (K/0.15)*length(cm));
            vals(vals < 1) = 1;
            vals(vals > 255) = 255; 
            col = cm( vals, :);
            col = [ permute(col, [3 1 2]); permute(col, [3 1 2])];
            z = zeros([ 1, size( coords,2)]);
            surface([coords(1,:);coords(1,:)],[coords(2,:);coords(2,:)],[z;z],col,...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',3);

        end
        % }}}
        
        % displayFeature3D {{{
        function ax = displayFeature3D( obj, ax, sizeZ, cm)
            if nargin < 4
                cm = hsv;
            end
             % Get (x,y,z) coordinates of the curve
            coords = obj.GetCoords();
            line( coords(1,:), coords(2,:), coords(3,:), obj.display{:} );
            x = [coords(1,:),NaN];
            y = [coords(2,:),NaN];
            z = [coords(3,:),NaN];
%             cm = hsv;
%             colormap(hsv)
            cols = round((z/sizeZ)*length(cm)); cols(end)=cols(end-1);
            cols( cols < 1) = 1; cols(cols>256) = 256;
%             cols = cols/256;
            col = 256*cm( cols, :);
            col = permute(col, [1 3 2]);
            
            p=patch(x,y,z,col/256,'FaceColor','Flat','EdgeColor','Flat','LineWidth',8);
        end
        % }}}

        % fillParams {{{
        function obj = fillParams( obj, sizeImage)
            obj.params = BundleNew.findVoxelsNearCurve( obj.GetCoords(), sizeImage, 6);
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

        % MakeModelFunction
        function ff = MakeModelFunction(obj)
            
            % Fourier model
            if contains(obj.model,'fourier')

                k = obj.curvature;
                a0 = k(1); k(1) = [];
                w = k(end); k(end) = [];
                aa = k(1:2:end);
                bb = k(2:2:end);
                ff = @(t) a0 + sum( aa.*cos( [1:length(aa)]*t*w)) + sum( bb.*sin( [1:length(bb)]*t*w));

            end

        end

        % GetCoords {{{
        function cc = GetCoords( obj, tmax, tmin)
            % get coordinates of curve from tmin to max. If not specified, get coords of entire curve
            if nargin < 2
                tmax = obj.L;
            end
            if nargin < 3
                tmin = 0;
            end
            if tmax > obj.L
                tmax = obj.L;
            end
            
            Rot = [ 0 -1; 1 0]; % Rotation matrix
            if obj.dim == 2
                tanVec = [ cos( obj.thetaInit(1)), sin( obj.thetaInit(1) )];
            elseif obj.dim == 3
                tanVec = [ cos( obj.thetaInit(1)), sin( obj.thetaInit(1) ), cos( obj.thetaInit(2) )];
            end
            % construct time vector
            t = linspace(0, tmax, ceil(1*tmax));  [~,idxStart] = min( abs(t-tmin));

            % Curvature function discretized in time
            acc = obj.MakeModelFunction();
            
            max_curvature = 0.05;
            
            % initialization
            xx = zeros( obj.dim, length(t) ); xx(:,1) = obj.origin';
            vv = zeros( numel(tanVec), length(t) ); vv(:,1) = tanVec;

            % iterate
            for jt = 2 : length(t)
                
                % Step size
                dt = t(jt)-t(jt-1);

                % New tangent : Get by adding normal vector to old tangent
                old_tan = vv(:,jt-1);
                acc_old = sign( acc(jt-1))*min([ abs(acc(jt-1)), max_curvature]);
                new_tan_XY = old_tan(1:2) + acc_old*dt*Rot*old_tan(1:2);
                vv(1:2,jt) = new_tan_XY/norm(new_tan_XY);
                
                if obj.dim == 3
                    vv(3,jt) = tanVec(3);
                end
                    
                xx(1:2,jt) = xx(1:2,jt-1) + dt*vv(1:2,jt);
                if obj.dim == 3
                    xx(3,jt) = xx(3,jt-1) + vv(3,jt)*dt;
                end
            end
            cc = xx(:,idxStart:end);

        end
        % }}}
        
        % GetMeanCurvature {{{
        function mean_K = GetMeanCurvature( obj)
            % get coordinates of curve from tmin to max. If not specified, get coords of entire curve
            tmax = obj.L;
            tmin = 0;
            
            % construct time vector
            t = linspace(0, tmax, ceil(10*tmax));

            % Acceleration function discretized in time
            % Curvature function discretized in time
            acc = obj.MakeModelFunction();
            mean_K = mean( abs(acc(t)) );

        end
        % }}}
        
        % GetMeanCurvature {{{
        function tip_K = GetTipCurvature( obj)
            % get coordinates of curve from tmin to max. If not specified, get coords of entire curve
            tmax = obj.L;
            tmin = 0;
            
            % construct time vector
            t = linspace(tmin, tmax, ceil(10*tmax));

            % Acceleration function discretized in time
            % Curvature function discretized in time
            acc = obj.MakeModelFunction();
            tip_K = acc(t(end));

        end
        % }}}
        % }}}

        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
            % Force features to lie inside the mask. This will
            % shorten any curves whose mtoc is within the mask.

            % Simulate image
            imFeat = obj.simulateFeature( size(mask) );
            
            % if escapes from mask, shorten until inside mask
            outside = Methods.CheckEscapeMask( imFeat, mask, 0.02);

            while outside

                % Shorten
                if obj.L-0.5 > obj.bounds.lb.L
                    obj.L = obj.L-0.5;
                else
                    break;
                end

                % Check again
                imFeat = obj.simulateFeature( size(mask) );
                outside = Methods.CheckEscapeMask( imFeat, mask);
                obj.GetLength();

            end
            
        end
        % }}}

        % GetStructInfo {{{
        function feat = GetStructInfo(obj)
            feat.type = type;
            feat.origin = obj.origin;
            feat.tanInit = obj.tanInit;
            feat.curvature = obj.curvature;
            feat.L = obj.L;
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

            S.model = obj.model;
            S.type = obj.type;
            S.dim = obj.dim;
            S.props2Fit = obj.props2Fit;
            S.amplitude = obj.amplitude;
            S.sigma = obj.sigma;
            S.origin = obj.origin;
            S.thetaInit = obj.thetaInit;
            S.curvature = obj.curvature;
            S.L = obj.L;
            S.display = obj.display;
            S.err_amplitude = obj.err_amplitude;
            S.err_sigma = obj.err_sigma;
            S.err_origin = obj.err_origin;
            S.err_thetaInit = obj.err_thetaInit;
            S.curvature = obj.curvature;
            S.err_L = obj.err_L;

        end
        % }}}
        
         % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
            % Get 2D projection of feature
            
            % Check object dimensionality
            if obj.dim == 2
                warning('object dimensionality is already 2')
            end
            
            obj.origin = obj.origin(1:2);
            obj.thetaInit = obj.thetaInit(1);
            obj.SetBounds();
            
        end
        % }}}
        
        % GetProjection3DSpecific {{{
        function obj = GetProjection3DSpecific( obj)
            % Get 3D projection of feature
            
            % Check object dimensionality
            if obj.dim == 3
                warning('object dimensionality is already 3')
            end
            
            obj.origin(3) = 7;
            obj.thetaInit(2) = 0;
            obj.SetBounds();
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            
            obj.origin(1:2) = obj2D.origin(1:2);
            obj.thetaInit(1) = obj2D.thetaInit(1);
            obj.sigma(1:2) = obj2D.sigma(1:2);
            obj.amplitude = obj2D.amplitude;
            obj.L = obj2D.L;
            obj.curvature = obj2D.curvature;
            
        end
        % }}}

        % SetBounds {{{
        function SetBounds( obj, pars)
           
            % origin
            tanVec = round( abs( 7*[cos(obj.thetaInit(1,1)), sin(obj.thetaInit(1,1)) ]));
            
            if obj.dim == 3
                ub.origin = [ obj.origin(1)+tanVec(1)+5, obj.origin(2)+tanVec(2)+5, obj.origin(3)-3];
                lb.origin = [ obj.origin(1)-tanVec(1)-5, obj.origin(2)-tanVec(2)-5, obj.origin(3)+3];
            elseif obj.dim == 2
                ub.origin = [ obj.origin(1)+tanVec(1)+5, obj.origin(2)+tanVec(2)+5];
                lb.origin = [ obj.origin(1)-tanVec(1)-5, obj.origin(2)-tanVec(2)-5];
            end
            
            % amplitude
            if nargin==2
                ub.amplitude = pars.amplitude.ub;
                lb.amplitude = pars.amplitude.lb;
                if obj.amplitude < lb.amplitude
                    obj.amplitude = 1.05*lb.amplitude;
                    disp('forcing curved MT amplitude above a lower threshold for Budding yeast')
                end
            else
                if ~isfield(obj.bounds, 'ub')
                    ub.amplitude = 1;
                    lb.amplitude = 0;
                else
                    ub.amplitude = obj.bounds.ub.amplitude;
                    lb.amplitude = obj.bounds.lb.amplitude;
                end
            end
            
            % sigma 
            if obj.dim == 3
                ub.sigma = [3.0 3.0 2.5];
                lb.sigma = [1.8 1.8 1.0];
            elseif obj.dim == 2
                ub.sigma = [3.0 3.0];
                lb.sigma = [1.8 1.8];
            end
            
            % L 
            ub.L = obj.L+20;
            lb.L = max( [6, obj.L-20]);

            % thetaInit
            if obj.dim == 3
                ub.thetaInit = obj.thetaInit + [0.3, 0.2];
                lb.thetaInit = obj.thetaInit - [0.3, 0.2];
            elseif obj.dim == 2
                ub.thetaInit = obj.thetaInit + [0.3];
                lb.thetaInit = obj.thetaInit - [0.3];
            end
            % normalVec
%             if length(obj.normalVec) == 1
%                 ub.normalVec = obj.normalVec + [0.02];
%                 lb.normalVec = obj.normalVec - [0.02];
%             elseif length(obj.normalVec) == 2
%                 ub.normalVec = obj.normalVec + [0.02, 0.0002];
%                 lb.normalVec = obj.normalVec - [0.02, 0.0002];
%             elseif length(obj.normalVec) == 3
%                 ub.normalVec = obj.normalVec + [0.02, 0.0002, 0.000002];
%                 lb.normalVec = obj.normalVec - [0.02, 0.0002, 0.000002];
%             end
            ub.curvature = obj.curvature + 0.02;
            lb.curvature = obj.curvature - 0.02;
            
            obj.bounds.lb = lb;
            obj.bounds.ub = ub;
            
        end
        % }}}
        
        % preOptimize {{{
        function obj = preOptimize(obj, imOrg, imBkg)
           
%             % optimize thetaInit and normalVec
%             res1 = []; 
%             deltaT = linspace( -0.5, 0.5, 10);  t0 = obj.thetaInit;
%             deltaN = linspace( -0.005, 0.005, 10); nv0 = obj.normalVec;
% %             figure; imagesc( max(imOrg,[],3)); axis equal; hold on;
%             for it = deltaT
%                 res2=[]; 
%                 for in = deltaN
%                     obj.thetaInit = t0 + [it, it]; 
%                     obj.normalVec = nv0 + [in, in];
%                     imSim = imBkg + obj.simulateFeature( size(imBkg));
%                     res2 = [ res2, sum( (imSim(:) - imOrg(:) ).^2 )];
%                     
% %                     cc = obj.GetCoords();
% %                     plot( cc(1,:), cc(2,:), 'LineWidth', 1);
%                     
%                 end
%                 res1 = [res1; res2];
%             end
%             [~,idT] = min( min(res1,[],2) ); obj.thetaInit = t0 + deltaT(idT); 
%             [~,idN] = min( min(res1,[],1) ); obj.normalVec = nv0 + deltaN(idN);
            

            % optimize amp
            if any( strcmp( obj.props2Fit, 'amplitude') )
                res = []; amps = linspace( 0, max( imOrg(:))-imBkg(1), 20); 
                A0 = obj.amplitude; 
                for ia = amps
                    obj.amplitude = ia; 
                    imSim = imBkg + obj.simulateFeature( size(imBkg));
                    res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
                end
                [~,idA] = min(res); 
                obj.amplitude = amps(idA);
            end

            % optimize sigma
%             if any( strcmp( obj.props2Fit, 'sigma') )
%                 % sx
%                 res = []; sx = linspace(2,3,10); sx0 = obj.sigma(1); sy0 = obj.sigma(2);
%                 for ix = sx
%                     obj.sigma(1) = ix; obj.sigma(2) = ix;
%                     imSim = imBkg + obj.simulateFeature( size(imBkg));
%                     res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
%                 end
%                 [~, idx] = min( res); 
%                 obj.sigma(1) = sx(idx);
%                 obj.sigma(2) = sx(idx);
%                 % sz
%                 if obj.dim == 3
%                     res = []; sz = linspace(1.0,2.5,10); sz0 = obj.sigma(3);
%                     for ix = sz
%                         obj.sigma(3) = ix;
%                         imSim = imBkg + obj.simulateFeature( size(imBkg));
%                         res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
%                     end
%                     [~, idx] = min( res); obj.sigma(3) = sz(idx);
%                 end
%             end
            
            % optimize lengths
%             if any( strcmp( obj.props2Fit, 'L') )
%                 % end 1
%                 res = []; l1 = linspace(obj.L,3*obj.L,10); l10 = obj.L;
%                 for ix = l1
%                     obj.L = ix;
%                     imSim = imBkg + obj.simulateFeature( size(imBkg));
%                     res = [ res, sum( (imSim(:) - imOrg(:) ).^2 )];
%                 end
%                 [~, idx] = min( res); 
%                 obj.L = l1(idx);
%             end
            %thr = multithresh(imOrg(:),2);
            %par.amplitude.lb = thr(1)*2;
            par.amplitude.lb = 1.2*median( imOrg(:) );
            par.amplitude.ub = max( imOrg(:));
            obj.SetBounds(par);
            
        end
        % }}}
    end

    methods ( Static = true )
        
        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'CurvedMT')
                error('incorrect type')
            end          
            
            obj = CurvedMT( 'fourier2',S.origin, S.thetaInit, S.curvature, S.L, S.amplitude, S.sigma, S.dim, S.props2Fit, S.display);
            try
                obj.err_origin = S.err_origin;
                obj.err_curvature = S.err_curvature;
                obj.err_thetaInit = S.err_thetaInit;
                obj.err_L = S.err_L;
                obj.err_amplitude = S.err_amplitude;
                obj.err_sigma = S.err_sigma;
            catch
                warning('errors not present')
            end
        end
        % }}}
        
        % findVoxelsNearCurve {{{
        function lineVox = findVoxelsNearCurve( coords, sizeImage, rad)
            
            dim = length( sizeImage);
            if size(coords,1) ~= dim
                error('Coord size doesnt match dim of image')
            end
            % Get pixels coordinates perp to direction of curve
            orient = atan2( coords(2,end)-coords(2,1) , coords(1,end)-coords(1,1) );
            nVec = [cos(orient+pi/2); sin(orient+pi/2) ];
            rads = -rad:rad;
            c2 = zeros( 2, size(coords,2)*length(rads));
            for jr = 1: length(rads)
                c2(:, 1+(jr-1)*size(coords,2):(jr)*size(coords,2)) = round( coords(1:2,:)+ rads(jr)*nVec);
            end
            % Add half circle near endpoints
            imt = zeros(sizeImage(1:2)); 
            imt(round(coords(2,1)),round(coords(1,1)))=1; imt(round(coords(2,end)),round(coords(1,end)))=1; 
            imt = imdilate( imt, strel('disk',2*rad)); idxEnd = find(imt);
            [ye, xe] = ind2sub( sizeImage(1:2), idxEnd);
            c2 = [c2 , [xe';ye']];
            pts = unique( c2', 'rows');
            
            % store voxel info
            if dim == 2
                lineVox.x = pts(:,1); lineVox.y = pts(:,2);
                lineVox.idx = sub2ind( sizeImage, pts(:,2), pts(:,1) );
            elseif dim == 3
                ps = zeros( sizeImage(3)*size(pts,1), 3);
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
        function coeffs = estimatePolyCoefficients( coords, order, t)
            % coords: D x N array (D: dimensionality, N: number of points
            % along curve)
            % order: This is the order of the polynomial: positive integer
            
            if length(order) ~= size(coords,1)
                error('dimensionality does not match')
            end

            % interpolate for better accuracy
            if nargin < 3
                t = linspace( 0, 1, size(coords,2));
            end
%             
%             xi = coords(1,:);
%             yi = coords(2,:);
            
            % fit polynomial of order polyOrder
%             if size(coords,1) == 2
%                 cc = [xi; yi];
%                 coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) };
%             elseif size(coords,1) == 3
%                 cc = [xi;yi;zi];
%                 zi = interp1( linspace(0,1, size(coords,2) ), coords(3,:), t, 'linear', 'extrap');
%                 coeffs = { polyfit( t, xi, order(1) ), polyfit( t, yi, order(2) ) , polyfit( t, zi, order(3) ) };
%             end
            
            coeffs = {}; cc=coords;
            for jc = 1 : size(cc,1)

                % fit using polyfit
                cf1 = polyfit(t, cc(jc,:), order(jc) );

                if order(jc) == 3
                    fitfunc = @(a,b,c,x) a.*(x.*x.*x) + b.*x.*x + c.*x + cc(jc,1);
                elseif order(jc) == 4
                    fitfunc = @(a,b,c,d,x) a.*(x.*x.*x.*x) + b.*x.*x.*x + c.*x.*x + d.*x + cc(jc,1);
                elseif order(jc) == 2
                    fitfunc = @(a,b,x) a.*(x.*x) + b.*x + cc(jc,1);
                elseif order(jc) == 1
                    fitfunc = @(a,x) a.*x + cc(jc,1);
                end
                ff = fit(t(:),cc(jc,:)', fitfunc, 'StartPoint',cf1(1:end-1) );
                cf2 = [ coeffvalues(ff), cc(jc,1) ];
                coeffs = {coeffs{:}, cf2};

            end
            
        end
        % }}}

    end

end

