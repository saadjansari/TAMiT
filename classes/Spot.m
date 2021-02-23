classdef Spot < BasicElement

    properties
        position
        err_position
    end

    methods ( Access = public )
       
        % Spot {{{
        function obj = Spot( position, amplitude, sigma, dim, props2Fit, display)
        % Spot : this is the constructor function for a Spot. This will be used for SPBs, MTOCs, and Kinetochores

            obj = obj@BasicElement( dim, amplitude, sigma, props2Fit, display, 'Spot');
            obj.position = position;
            obj.SetBounds();

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels, ub,lb] = getVec( obj, props2get)

            % sample props2get
            if nargin==1
                props2get = {'position', 'amplitude', 'sigma'};
            end
            
            % make sure props input matches properties defined in class
            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
            end

            % Get vector of Properties
            if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                temps.position = obj.position(3);
                temps.sigma = obj.sigma(3);
                temps.amplitude = obj.amplitude;
                temps.bounds.ub.position = obj.bounds.ub.position(3);
                temps.bounds.ub.sigma = obj.bounds.ub.sigma(3);
                temps.bounds.ub.amplitude = obj.bounds.ub.amplitude;
                temps.bounds.lb.position = obj.bounds.lb.position(3);
                temps.bounds.lb.sigma = obj.bounds.lb.sigma(3);
                temps.bounds.lb.amplitude = obj.bounds.lb.amplitude;
            end
            vec = []; ub = []; lb = [];
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
            if isempty(ub)
                clearvars ub lb
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels, errBoolean)

            if nargin < 4
                errBoolean = 0;
            end
            
            props2find = {'position', 'amplitude', 'sigma'};

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
                    continue
                end
                
                if obj.dim == 3 && strcmp(obj.fit, 'zonly')
                    obj.(propCurr) = obj.( props2find{ jProp});
                    obj.( propCurr )(1+end-length(idxProp):end) = vec( idxProp);
                    
                else
                    if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                        error( 'absorbVec: length of vector props to absorb does not match the old property size')
                    end

                    % Set final property
                    obj.( propCurr ) = vec( idxProp);
                end

            end

        end
        % }}}
        
        % simulateFeature {{{
        function imageOut = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: must provide size of box to simulate feature in')
            end

            % Simulate a gaussian spot
            imageOut = zeros( sizeImage);
            if isfield( obj.params, 'idx')
                if obj.dim == 2
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Spot2', 'Pos', obj.position, 'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y);
                elseif obj.dim == 3
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Spot3', 'Pos', obj.position, 'Idx', obj.params.idx, 'X', obj.params.x, 'Y', obj.params.y, 'Z', obj.params.z);
                end
                imageOut = obj.amplitude * mat2gray( imGraph);
            else
                if obj.dim == 2
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Spot2', 'Pos', obj.position);
                elseif obj.dim == 3
                    [imGraph, ec, err] = DrawGaussian( obj.sigma, imageOut, 'Spot3', 'Pos', obj.position);
                end
                imageOut = obj.amplitude * mat2gray( imGraph);
            end

            % What to do if there is an error in Z?
            % Options:
            % 1. Increase intensity at closest 3D point
            if err > 0
                errorPlane = imageOut(:,:,ec);
                idx = find( errorPlane == max( errorPlane(:) ) );
                [yidx,xidx] = ind2sub( size( errorPlane), idx);
                imageOut( yidx,xidx, ec) = imageOut( yidx,xidx, ec)*(1 +err);
            end
            
%             imageFeat = imageOut;
            obj.imageSim = imageOut;
%             imageOut( imageOut < imageFeat) = imageFeat( imageOut < imageFeat);
%             imageOut = imageFeat + imageOut;

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, sizeZ)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.image, [], 3) ); colormap gray; axis equal;
            end
            
            if obj.dim==3 && nargin==3
                
                cm = hsv;
                col = cm( round((obj.position(3)/sizeZ)*length(cm)), :);
                h1 = plot( obj.position(1), obj.position(2), 'Marker', 'o', 'MarkerSize',10, 'Color', 'w', 'LineWidth', 2);
                set(h1, 'markerfacecolor', col);
                % colorbar('Ticks',linspace(0,1,sizeZ),'TickLabels',1:sizeZ)
                
            else
                % Create the spot to display
                plot( obj.position(1), obj.position(2), obj.display{:} );
            end

        end
        % }}}
        
        % displayFeature3D {{{
        function ax = displayFeature3D( obj, ax,sizeZ)
            cm = hsv;
            col = cm( round((obj.position(3)/sizeZ)*length(cm)), :);
            line( obj.position(1), obj.position(2), obj.position(3), 'Color', col, ...
                'Marker','o', 'MarkerSize',20,'LineWidth',10,'MarkerFaceColor',col );
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

            plot( obj.position(3), obj.position(1), obj.display{:} ) 

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

            plot( obj.position(3), obj.position(2), obj.display{:} ) 

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
        function obj = fillParams( obj, sizeImage)
            obj.params = Spot.findVoxelsNearSpot( obj.position, sizeImage, 10);
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
            S.err_position = obj.err_position;
            S.err_amplitude = obj.err_amplitude;
            S.err_sigma = obj.err_sigma;
        end
        % }}}
        
        % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
            % Get 2D projection of feature
            
            % Check object dimensionality
            if obj.dim == 2
                warning('object dimensionality is already 2')
            end
            
            obj.position = obj.position(1:2);
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
            
            obj.position(3) = 1;
            obj.SetBounds();
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            
            obj.position(1:2) = obj2D.position(1:2);
            obj.sigma(1:2) = obj2D.sigma(1:2);
            obj.amplitude = obj2D.amplitude;
            
        end
        % }}}
        
        % GetStructInfo {{{
        function feat = GetStructInfo(obj)
            feat.type = obj.type;
            feat.position = obj.position;
            feat.amplitude = obj.amplitude;
            feat.sigma = obj.sigma;
            feat.ID = obj.ID;
        end
        % }}}
        
        % SetBounds {{{
        function SetBounds( obj)
           
            % amplitude
            ub.amplitude = 1;
            lb.amplitude = 0;
            
            % sigma % positions
            if obj.dim == 3
                if strcmpi(obj.label, 'spb')
                	ub.sigma = [3.0 3.0 2.0];
                else
                    ub.sigma = [2.0 2.0 2.0];
                end
                lb.sigma = [1.2 1.2 1.0];
                ub.position = [150 150 7];
                lb.position = [1 1 1];
            elseif obj.dim == 2
                if strcmpi(obj.label, 'spb')
                	ub.sigma = [3.0 3.0];
                else
                    ub.sigma = [2.0 2.0];
                end
                lb.sigma = [1.2 1.2];
                ub.position = [150 150];
                lb.position = [1 1];
            end
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
            obj.SetBounds();
            
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
            try
                obj.err_position = S.err_position;
                obj.err_amplitude = S.err_amplitude;
                obj.err_sigma = S.err_sigma;
            catch
                warning('errors not present')
            end
        end
        % }}}
        
        % findVoxelsNearSpot {{{
        function spotVox = findVoxelsNearSpot( pos, sizeImage, radDilate)
            
            dim = numel( pos );
            if dim ~= 2 && dim ~= 3
                error('findVoxelsNearSpot: dim must be 2 or 3')
            end
            
            % Draw the spot in empty space
            imSpot = zeros( sizeImage);
            if dim == 2
                imSpot( round(pos(2)), round( pos(1) ) ) = 1;
            elseif dim == 3
                imSpot( round(pos(2)), round( pos(1) ), round( pos(3) ) ) = 1;
            end
            
            % dilate the line with a sphere of large size
%             imSpot = imdilate( imSpot, strel( 'sphere', radDilate) );
            imSpot = imdilate( max( imSpot, [], 3), strel('disk', radDilate) );
            if dim == 3
                imSpot = repmat( imSpot, 1, 1, sizeImage(3) );
            end
            
            % return voxel indices for the dilated line
            spotVox.idx = find( imSpot(:) );
            if dim == 2
                [spotVox.y, spotVox.x] = ind2sub( size(imSpot), spotVox.idx);
            elseif dim == 3
                [spotVox.y, spotVox.x, spotVox.z] = ind2sub( size(imSpot), spotVox.idx);
            end
            
            
        end
        % }}}
        
    end
end
