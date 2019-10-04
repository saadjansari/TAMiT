classdef Cell < handle & matlab.mixin.Copyable
% This is a generic cell which can be specialized
    properties
        type % type of cell (e.g. interphase or mitotis). 
        imageData % (X,Y,Z,T,C)
        species
        strain % any specific strain information
        featuresInChannels
        channelsToFit
        featureList % cell array containing features. this is of size numberOfChannels x Time.
        params
        featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    end

    methods ( Access = public )
        
        % Cell {{{
        function obj = Cell( imageData, featuresInChannels, channelsToFit, params, varargin)
        % Cell : constructor function for Cell superclass
       
            % Store essentials
            obj.imageData = imageData; % imageData object
            obj.featuresInChannels = featuresInChannels;
            obj.channelsToFit = channelsToFit;
            obj.params = params;

            % parse options
            opts = parseArgs( varargin{:} );
            obj.type = opts.Type;
            obj.species = opts.Species;
            obj.strain = opts.Strain;

            % Initialize the feature list
            obj.featureList = cell( length( obj.channelsToFit), size( imageData.GetImage, 4) );

            % Save Directory creation 
            if exist( obj.params.saveDirectory) ~= 7
                mkdir( obj.params.saveDirectory)
            end
            
            % Define a global counter for objects
            global counter; counter = -1;

            Cell.PrintInfo( obj);

            % parseArgs {{{
            function opts = parseArgs( varargin)
               
                % default 
                defaultType = 'generic';
                defaultSpecies = 'Pombe';
                defaultStrain = 'WildType';

                % Input Parser
                p = inputParser;
                addParameter( p, 'Type', defaultType);
                addParameter( p, 'Species', defaultSpecies);
                addParameter( p, 'Strain', defaultStrain);

                parse( p, varargin{:});
                opts = p.Results;

            end
            % }}}
            
        end
        % }}}
        
        % FindFeatures {{{
        function obj = FindFeatures( obj)

            for jChannel = 1 : length(obj.channelsToFit)

                cChannel = obj.channelsToFit( jChannel);
                disp( upper( sprintf( ['Channel %d = ' obj.featuresInChannels{jChannel}], cChannel ) ) )
                disp('-----------------------------------------------------------------------')

                obj.FindFeaturesChannel( jChannel);

            end

        end
        % }}}

        % FindFeaturesChannel {{{
        function obj = FindFeaturesChannel( obj, jChannel)

            % get actual channel number based on channel index
            cChannel = obj.channelsToFit( jChannel);

            % Get lifetime for this cell
            lifetime = obj.imageData.GetLifetime;

            for jTime = lifetime(1) : lifetime(2)

                disp( sprintf( 'C-%d Time = %d', jChannel, jTime) )

                % Initialize for this CT step
                params = obj.params;
                params.channelTrue = cChannel;
                params.channelIdx = jChannel;
                params.time = jTime;
                params.saveDirectory = [ obj.params.saveDirectory , filesep, sprintf( 'C%d_T%d', cChannel, jTime) ];
                mkdir( params.saveDirectory);

                % Find features for this CT frame
                obj.FindFeaturesFrame( params) 

            end
            disp('-----------------------------------------------------------------------')
            
        end
        % }}}
        
        % FindFeaturesFrame {{{
        function FindFeaturesFrame( obj, parameters)
            % Find and fit features for a single CT frame
            
            % Get the image for this frame
            Image = obj.imageData.GetImage();
            Image2Fit = Image(:,:,:, parameters.time, parameters.channelTrue);

            % check if frame is bad. In this case, we'll use the previous fitted frame and then skip fitting
            if Cell.checkFrameViability( Image2Fit, parameters.time)~= 0
                obj.featureList{ parameters.channelIdx, parameters.time} = obj.featureList{ parameters.channelIdx, parameters.time-1};
                fprintf( 'Encountered a bad frame. Using previous frame features and skipping fitting...\n')

                % save Final fit
                p.channel = parameters.channelTrue; p.time = parameters.time; p.fitScope = 'global';
                Cell.saveFinalFit( Image2Fit, obj.featureList{ parameters.channelIdx, parameters.time}, p);
                return
            end

            % Estimate the features based on an estimation routine (defined in specialized sub-class )
            % Good estimation is key to good optimization in low SnR images
            disp('Estimating features...')
            obj.EstimateFeatures( parameters.time, parameters.channelTrue, parameters.channelIdx); 
            mainFeature = obj.featureList{ parameters.channelIdx , parameters.time};
            obj.syncFeatureMap( parameters.channelIdx, parameters.time);

            % Prepare fit params
            params = obj.params.fit;
            params.channel = parameters.channelTrue;
            params.time = parameters.time;
            params.saveDirectory = parameters.saveDirectory;

            % Fit Features via Fit Engine.
            fitEngine = FitEngine( Image2Fit, mainFeature, params);
            fitEngine.Optimize();

            obj.syncFeatureMap( parameters.channelIdx, parameters.time);
            close all

        end
        % }}}
        
        % simulateCell {{{
        function imageOut = simulateCell( obj, cChannel, cTime)
            % Simulates XYZ image of a cell with a mainFeature
            
            mainFeature = obj.featureList{cChannel, cTime};

            imageOut = 0*mainFeature.image;

            % simulate the main feature 
            imageOut = mainFeature.simulateFeature( imageOut);

        end
        % }}}
        
        % syncFeatureMap {{{
        function obj = syncFeatureMap( obj, channel, time)
            % Updates the counter and assigns IDs to all subfeatures
           
            global COUNTER 
            mapGlobal = obj.featureMap;
            keysGlobal = keys( mapGlobal);
            valsGlobal = values( mapGlobal);
            mainFeature = obj.featureList{ channel, time}; 

            % Remove any entries for this channel and time
            for jKey = 1 : length( keysGlobal)
                if all( valsGlobal{ jKey}(1:2) == [channel time])
                    remove( mapGlobal, keysGlobal( jKey) );
                end
            end

            % Get current counter value
            COUNTER = uint32(max( cell2mat( keysGlobal ) ) + 1);
            if isempty( COUNTER), COUNTER = uint32(1); end
            
            % Assign an ID to its main feature. and prompt it to assign IDs to its sub features
            mainFeature.ID = COUNTER; 
            COUNTER = COUNTER + 1;
            mainFeature.syncFeatures();

            % Update global map
            % Get mainFeature local map, append current channel and time and add it to our map
            mapLocal = mainFeature.featureMap;
            keysLocal = cell2mat( keys( mapLocal) );
            valsLocal = values( mapLocal);

            for jKey = 1 : length( keysLocal)
                mapGlobal( keysLocal( jKey) ) = [ channel time valsLocal{ jKey}];
            end
            obj.featureMap = mapGlobal;
            mainFeature.featureMap;
            
        end
        % }}} 

    end
    
    methods ( Static = true , Access = public )
        
       % Common Functions {{{
       
        % img2double {{{
        function imgOut = img2double( imgIn)

            mx = double( max( imgIn(:) ) );
            mn = double( min( imgIn(:) ) );
            dif = mx-mn;
            imgOut = mn + dif*im2double( imgIn);

        end
        % }}}
        % img2uint16 {{{
        function imgOut = img2uint16( imgIn)

            mx = double( max( imgIn(:) ) );
            mn = double( min( imgIn(:) ) );
            dif = mx-mn;
            imgOut = mn + dif*im2double( imgIn);

        end
        % }}}
        % findAmplitudeAlongLine {{{
        function amplitude = findAmplitudeAlongLine( imageIn, startPos, endPos)

            dim = length( size( imageIn) );
            if dim ~= length( startPos) || dim ~= length( endPos)
                error( 'findAmplitudeAlongLine: dimensions of input arguments does not match') 
            end

            imageIn = im2double( imageIn);

            LineLength = norm( abs(endPos - startPos) );

            numPoints = round( LineLength*5);
            if dim == 2
                xq = linspace( startPos(1), endPos(1), numPoints);
                yq = linspace( startPos(2), endPos(2), numPoints);
                amplitude = interp2( imageIn, yq, xq, 'linear');

            elseif dim == 3
                xq = linspace( startPos(1), endPos(1), numPoints);
                yq = linspace( startPos(2), endPos(2), numPoints);
                zq = linspace( startPos(3), endPos(3), numPoints);
                amplitude = interp3( imageIn, xq, yq, zq, 'linear');

            else
                error('findAmplitudeAlongLine: dimensionality of data must be 2D or 3D')
            end

        end
        % }}}
        % radIntegrate2D {{{
        function [phiIntensity, phiValues] = radIntegrate2D( imageIn, startPoint)

            if length( size(imageIn) ) ~= 2
                error( 'radIntegrate2D : input image must be 2-dimensional')
            end

            % Radially integrate around startPoint 
            PhiStep = deg2rad(2);
            rmin = 1;
            rmax = 10;
            numVoxelsX = size( imageIn, 2);
            numVoxelsY = size( imageIn, 1);
            x0 = startPoint(1);
            y0 = startPoint(2);

            % Define the step sizes to use for the angular sweep in phi. 
            PhiStep = deg2rad(2);
            rmin = 1;
            rmax = 10;
            phiValues = 0 : PhiStep : 2*pi - PhiStep;

            % Pre-allocate array for storing intensity data with varying phi.
            phiIntensity = zeros( 1, length( phiValues) );

            method = 'interp2';
            switch method

            % method: interp2 {{{
            case 'interp2'
                for jPhi = 1 : length( phiValues)

                    % Find the coords along this line
                    X1 = x0 + ( [rmin:0.5:rmax] .* cos( phiValues( jPhi) ) );
                    Y1 = y0 + ( [rmin:0.5:rmax] .* sin( phiValues( jPhi) ) ); 

                    % Identify coords that exceed the size of the image, and remove them
                    rmIdx = union( find( X1 < 1 | X1 > numVoxelsX) , find( Y1 < 1 | Y1 > numVoxelsY) );
                    X1( rmIdx) = [];
                    Y1( rmIdx) = [];
                    
                    % line Intensity along the line
                    lineInt = interp2( im2double( imageIn), X1, Y1, 'linear');
                    phiIntensity( jPhi) = mean( lineInt);
                    
                end
            % }}}

            % method: convexhull {{{
            case 'convexhull'
                for jPhi = 1 : length( phiValues)
                    % Extract the relevant angles for drawing a non-zero angle
                    phi1 = phiValues( jPhi) - PhiStep;
                    phi2 = phiValues( jPhi) + PhiStep;

                    % Find the coordinates of the extreme points on the arc length
                    % of this pizza slice.
                    X1 = x0 + ( rmax * cos( [ phi1, phi2]) );
                    Y1 = y0 + ( rmax * sin( [ phi1, phi2]) ); % -ve because matlab orientation is reversed( origin is at top left)

                    % find minimum points at distance rmin away
                    X2 = x0 + ( rmin * cos( [ phi1, phi2]) );
                    Y2 = y0 + ( rmin * sin( [ phi1, phi2]) ); 

                    X = [ X1, X2];
                    Y = [ Y1, Y2];

                    xRd = []; yRd = [];
                    for jPt = 1 : length(X)
                        [~, xidx] = min( abs( xVecNew - round( X(jPt), 1) ) );
                        [~, yidx] = min( abs( xVecNew - round( Y(jPt), 1) ) );
                        xRd = [ xRd, xidx ];
                        yRd = [ yRd, yidx ];
                    end

                    % Initialize the mask. We'll turn on the (X,Y) pixels and draw a convex
                    % hull around it.
                    imMask = zeros( size( imSub) );
                    imMask( sub2ind( size( imSub), yRd, xRd ) ) = 1;

                    % Draw a convex hull and uncover the pixel values of interest
                    imMask = bwconvhull( imMask);

                    imMasked = imMask .* imSub;
            %         imshowpair(imSub, imMask); pause(0.25)
                    
                    % Sum up values and store them
                    phiIntensity( jPhi) = sum( imMasked(:) ) / sum(imMask(:));
                    
                end
                % }}}

            end

            % We smooth the angular intensity to enable easy peak finding
            phiIntensity = imgaussfilt( phiIntensity, 3);
           
            % find idx of minimum intensity and shift the array so no peaks are on the edge
            [~, idxShift] = min( phiIntensity);
            phiIntensity = circshift( phiIntensity, -idxShift);
            phiValues( 1:idxShift) = phiValues(1:idxShift) + 2*pi;
            phiValues = circshift( phiValues, -idxShift);
%             figure; plot( phiValues, phiIntensity, 'r-')

        end
        % }}}
        % phiIntegrate2D {{{
        function [radIntensity, radValues] = phiIntegrate2D( imageIn, startPoint)

            if length( size(imageIn) ) ~= 2
                error( 'radIntegrate2D : input image must be 2-dimensional')
            end

            % Radially integrate around startPoint 
            radStep = 0.25;
            rmin = 1;
            rmax = 150;
            radValues = rmin: radStep : rmax;

            numVoxelsX = size( imageIn, 2);
            numVoxelsY = size( imageIn, 1);
            x0 = startPoint(1);
            y0 = startPoint(2);

            % Define the step sizes to use for the angular sweep in phi. 
            PhiStep = deg2rad(2);
            phiValues = 0 : PhiStep : 2*pi - PhiStep;

            % Pre-allocate array for storing intensity data with varying phi.
            phiIntensity = zeros( 1, length( phiValues) );

            % Get coordinates of circle perimeter for varying radia 
            Xcirc = x0 + rValues'.*cos( phiValues);
            Ycirc = y0 + rValues'.*sin( phiValues);

            % Get intensity values of these scattered points in 2D space
            Intcirc = interp2( im2double( imageIn), Xcirc, Ycirc, 'linear', 0);
           
            % Sum over phi to get a vector of Intensity values at varying radia
            radIntensity = sum( Intcirc, 2)';

        end
        % }}}
        % drawGaussianPoint3D {{{
        function imagePoint = drawGaussianPoint3D( pos, sigma, imageIn, idx, x, y, z)
            % Draws a gaussian point using an analytical framework

            % Check that everything is in 3D
            if numel( size(imageIn) ) ~= 3 || length(pos)~=3 || length(sigma)~=3
                error('drawGaussianLine3D: all data must be 3-dimensional')
            end

            if nargin < 4
                idx = 1 : numel( imageIn);
            end
            if nargin < 7
                [y, x, z] = ind2sub( size( imageIn), idx);
            end

            x0 = pos(1); y0 = pos(2); z0 = pos(3);
            sx = sigma(1); sy = sigma(2); sz = sigma(3);

            % Amplitudes:
            ExpX = exp( -0.5*( (x-x0)./sx).^2 );
            ExpY = exp( -0.5*( (y-y0)./sy).^2 );
            ExpZ = exp( -0.5*( (z-z0)./sz).^2 );
            IntValues = ExpX .* ExpY .* ExpZ;
           
            functionCheck = 0;
            if functionCheck
                disp( sprintf('Total Voxels = %d', numel( IntValues) ) )
                disp( sprintf('Nan Voxels = %d', sum( isnan( IntValues(:) ) ) ) )
                disp( sprintf('Inf Voxels = %d', sum( isinf( IntValues(:) ) ) ) )
            end
            IntValues( isnan( IntValues) ) = min( IntValues(:) );
            IntValues( IntValues == Inf) = min( IntValues(:) );
            
            % Initiliaze the volume and set the appropriate indices to these values
            imagePoint = 0*imageIn;
            imagePoint( idx) = IntValues;

            if functionCheck
                dispImg( imagePoint);
            end
                    
        end
        % }}}
        % drawGaussianLine3D {{{
        function imageLine = drawGaussianLine3D( startPos, endPos, sigma, imageIn, idx, x, y, z)
            % Draws a gaussian straight line using an analytical framework

            % Check that everything is in 3D
            if numel( size(imageIn) ) ~= 3 || length(startPos)~=3 || length(endPos)~=3 || length(sigma)~=3
                startPos
                endPos
                sigma
                size(imageIn)
                error('drawGaussianLine3D: all data must be 3-dimensional')
            end

            if nargin < 5
                idx = 1 : numel( imageIn);
            end
            if nargin < 8
                [y, x, z] = ind2sub( size( imageIn), idx);
            end

            x0 = startPos(1); y0 = startPos(2); z0 = startPos(3);
            x1 = endPos(1); y1 = endPos(2); z1 = endPos(3);
            sx = sigma(1); sy = sigma(2); sz = sigma(3);

            % First lets parameterize this line with a parameter t. We'll find the
            % slopes of the line in each spatial dimension along with any offset
            % At t = 0, x = x0, y = y0, z = z0
            % At t = 1, x = x1, y = y1, z = z1
            t0 = 0;
            t1 = 1;
            
            % find slopes
            mx = (x1 - x0) / ( t1 - t0);
            my = (y1 - y0) / ( t1 - t0);
            mz = (z1 - z0) / ( t1 - t0);
            
            % find offsets, these are just the start points at t = 0;
            cx = x0;
            cy = y0;
            cz = z0;
            
            % We will integrate a 3D gaussian with respect to the parameter t going
            % from 0 to 1.
            
            % We've done the analytical integration in mathematica and will use the
            % result here:
            
            % Amplitudes:
            AmpExp = - ( 1 / sqrt(mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2 ) );
            
            AmpErf = sqrt( pi/2 ) * sx * sy * sz;
            
            % Exponential factors:
            ExpDenom = 2 * ( mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2);
            
            Exp1 = ( - ( cx^2 * mz^2 * sy^2 + cz^2 * (my^2 * sx^2 + mx^2 * sy^2) + ...
                cx^2 * my^2 * sz^2 + cy^2 * (mz^2 * sx^2 + mx^2 * sz^2) ) / ExpDenom );
            
            Exp2 = ( - ( - 2 * cx * mz^2 *sy^2 * x - 2 * cx * my^2 * sz^2 * x + ...
                mz^2 * sy^2 * x.^2 + my^2 * sz^2 * x.^2 + 2 * cx * mx * my * sz^2 * y - ...
                2 * mx * my * sz^2 * x .* y + mz^2 * sx^2 * y.^2 + mx^2 * sz^2 * y.^2 ) / ExpDenom  );
            
            Exp3 = ( - ( 2 * cx * mx * mz * sy^2 * z - 2 * mx * mz * sy^2 * x .* z - ...
                2 * my * mz * sx^2 * y .* z + my^2 * sx^2 * z.^2 + mx^2 * sy^2 * z.^2  ) / ExpDenom );
            
            Exp4 = ( - ( - 2 * cz * ( cy * my * mz * sx^2 + cx * mx * mz * sy^2 - ...
                mx * mz * sy^2 * x - my * mz * sx^2 * y + my^2 * sx^2 * z + mx^2 * sy^2 * z )  ) / ExpDenom );
            
            Exp5 = ( - ( - 2 * cy * (cx *mx *my *sz^2 - mx * my * sz^2 * x + ...
                mx^2 * sz^2 * y + mz * sx^2 * (mz * y - my * z)) ) / ExpDenom );
           
            SumExp = Exp1 + Exp2 + Exp3 + Exp4 + Exp5;
            
        %     Exp5( Exp5 == Inf) = realmax;
            
            % Erf factors;
            ErfDenom = ( sqrt(2) * sx * sy * sz * sqrt( mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2 ) );
            
            Erf1 = erf( ( cz * mz * sx^2 * sy^2 + cy * my * sx^2 * sz^2 + cx * mx * sy^2 * sz^2 - ...
                mx * sy^2 * sz^2 * x - my * sx^2 * sz^2 * y - mz * sx^2 * sy^2 * z ) / ErfDenom );
            
            Erf2 = erf( ( cz * mz * sx^2 * sy^2 + mz^2 * sx^2 * sy^2 + ...
                sz^2 * (cy * my * sx^2 + my^2 * sx^2 + mx * sy^2 * (cx + mx - x) - my * sx^2 * y) - ...
                mz * sx^2 * sy^2 * z ) / ErfDenom );

            % Find the intensity values for these provided query x,y,z coordinates
            IntValues = AmpExp .*exp( SumExp) .* ...
                AmpErf .* ( Erf1 - Erf2 );
           
            functionCheck = 0;
            if functionCheck
                disp( sprintf('Total Voxels = %d', numel( IntValues) ) )
                disp( sprintf('Nan Voxels = %d', sum( isnan( IntValues(:) ) ) ) )
                disp( sprintf('Inf Voxels = %d', sum( isinf( IntValues(:) ) ) ) )
                figure; imagesc( max( imageLine, [], 3) );
                figure; imagesc( max( imageLine, [], 3) );
            end

            IntValues( isnan( IntValues) ) = min( IntValues(:) );
            IntValues( IntValues == Inf) = min( IntValues(:) );
            
            % Initiliaze the volume and set the appropriate indices to these values
            imageLine = 0 * imageIn;
            imageLine( idx) = IntValues;

        end
        % }}}
            % checkFrameViability {{{
            function status = checkFrameViability( Image2Fit, jTime);

                status = 0;
                % Check if frame is all super-bright( i.e. voxel variance is close to zero)
                maxVox = max( Image2Fit(:) );
                minVox = min( Image2Fit(:) );

                if maxVox == minVox
                    status = 1;
                    disp( sprintf( 'Unviable frame %d, Skipping it...', jTime) )
                end

            end
            % }}}
            
       % }}}
       
        % Feature Estimation {{{
        % find_manyLines {{{
        function [Lines, ax] = find_manyLines( imageIn, startPoint, ax)

            imageIn = im2double( imageIn);
            image2D = max( imageIn, [], 3);

            if nargin < 3 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end

            % radially integrate in phi
            [ phiIntensity, phiValues] = Cell.radIntegrate2D( image2D, startPoint);
%             figure; plot( phiValues, phiIntensity, 'r-');

            % find peaks in Phi Intensity
            imVals = image2D( image2D ~= 0);
            mtBkg = median( imVals );
            minPkHeight = mtBkg + std( imVals );
            
            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, 'MinPeakHeight', minPkHeight );
            peakPhi = mod( peakPhi, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )

            Lines = {};
            % Find length of lines
            for jLine = 1 : length( peakPhi)

                if nargout > 1
                    [ cLine, ax] = Cell.find_singleLine3D( imageIn, startPoint, peakPhi( jLine), ax );
                else
                    cLine = Cell.find_singleLine3D( imageIn, startPoint, peakPhi( jLine) );
                end

                if ~isempty( cLine)
                    Lines = { Lines{:}, cLine}; 
                end

            end
%             error('stop here')
        end
        % }}}
        % find_singleLine3D {{{
        function [Line, ax] = find_singleLine3D( imageIn, startPoint, startAngle, ax)
            % This works with both 3D and a 2D input input. The outputs vary for 2D and 3d inputs

            if nargin < 4 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( image2D, [], 3) ); colormap gray; axis equal; hold on;
            end

            dim = length( size( imageIn) );
            [image2D, idx3Max] = max( imageIn, [], 3);
            
            startPoint2D = startPoint(1:2);

            if nargout > 1
                [Line, ax] = Cell.find_singleLine2D( image2D, startPoint2D, startAngle, ax);
            else
                Line = Cell.find_singleLine2D( image2D, startPoint2D, startAngle);
            end

            if isempty( Line)
                return
            end

            if dim == 3
                % Find Z position for this Line, find theta angle for this Line

                % Approach : from the computed MIP, we also have access to the respective
                % z-indices that correspond to each maximum value. For each guessed
                % microtubule, we'll point it in its direction and get the local indices to
                % the nearest pixel. With these indices we can compute what theta must be.

                % define start radius for theta determination
                rmin = 1;

                % Define the radia up to which you'll check the z-positions
                RadVec = rmin : Line.length;
                ThetaVec = 0*RadVec;

                % for every possible point on the Line before it ends
                % final coords
                x1 = startPoint(1) + RadVec * cos( startAngle );
                y1 = startPoint(2) + RadVec * sin( startAngle );

                % get the z coord of the maxidx at this position, and calculate theta based on this value.
                z0 = startPoint(3);
                z1 = interp2( idx3Max, x1, y1, 'Linear');
                zdist = abs( z1 - z0);
                theta = atan( abs( z1 - z0) ./ RadVec);
               
                % Based on how theta is defined, starting at the positive z-axis
                % and sweeping pi radians until it hits the negative z-axis, we'll
                % find theta depending on the location of z1 above or below the SPB
                % in the z-dimension
                theta( z1 >= z0) = pi/2 - theta( z1 >= z0);
                theta( z1 < z0) = pi/2 + theta( z1 < z0);
                Line.theta = mean(theta);
                
                z1 = z0 + Line.length*cos( Line.theta);
                z1( z1 > size(imageIn, 3) ) = size( imageIn, 3);
                z1( z1 < 1) = 1;
                Line.startPosition = [ Line.startPosition , z0 ];
                Line.endPosition = [ Line.endPosition , z1 ];
            end

        end
        % }}}
        % find_singleLine2D {{{
        function [Line, ax] = find_singleLine2D( image2D, startPoint, startAngle, ax)

            % Input Checks
            if length( startPoint) ~= 2
                error('find_singleLine2D: input argument startPointmust be of length 2')
            end

            if nargin < 4 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( image2D, [], 3) ); colormap gray; axis equal; hold on;
            end

            Lmin = 5;

            % This will look at radial amplitude alon the Line to find when the intensity falls below a threshold
            threshIntensity = median( image2D( image2D ~= 0) ) + 1.5*std( image2D( image2D ~= 0) );
%             median( image2D( image2D ~= 0) ) 
%             std( image2D( image2D ~= 0) )
%             error('stop here')

            % Define what radia to use:
            rmin = 1;
            rmax = 20;
            radialValues = rmin : 0.5 : rmax;
            radialIntensity = 0* radialValues;
            x0 = startPoint( 1);
            y0 = startPoint( 2);
                
            for jRad = 1 : length( radialValues)

                jLen = radialValues( jRad);
                % Extract the relevant angles for drawing a non-zero angle

                % Find the coordinates of the extreme points on the arc length of this.
                X = x0 + ( jLen * cos( startAngle) );
                Y = y0 + ( jLen * sin( startAngle) );

                % Get the intensity value at this position
                radialIntensity( jRad) = interp2( image2D, X, Y, 'Linear');
                
            end

            % get the length value of this Line by finding the length along
            % the Line when the intensity drops below a certain threshold.
%             radialIntensity = smooth( radialIntensity, 3);

%             figure; plot(radialValues, radialIntensity); 
        
            [ ~, IdMax] = find( radialIntensity < threshIntensity, 1, 'first');
            if all(radialIntensity > threshIntensity), IdMax = length(radialIntensity); end

            if isempty( IdMax) || radialValues( IdMax) < Lmin, 
                Line = [];
                return
            else 
                % Store the corresponding properties for the Line 
                Line.length = radialValues( IdMax);
                Line.phi = startAngle;
                Line.startPosition = startPoint;
                Line.endPosition = startPoint + Line.length * [cos(startAngle), sin(startAngle)];
                if nargout > 1
                    line( [Line.startPosition(1), Line.endPosition(1)], [Line.startPosition(2), Line.endPosition(2)], 'Color', 'r', 'LineWidth', 3); 
                end

            end

        end
        % }}}
        % }}}
            
        % Graphics {{{
        
        % plot_midFit {{{
        function stop = plot_midFit( vec, optim, state, fitInfo )
            % Plots during Lsqnonlin fitting
            % Inputs :  
            %   OPTIMVALUES: properties are funccount, iteration
            %   STATE: string representing current state of optimization engine. Values are init, iter, done 
            %   STOP: A boolean to stop the algorithm.
            %   All other inputs available to error function are also available. Currently: ImageOrg( original image), fitInfo and parameters

            stop = false;

            ImageOrg = fitInfo.image;
            nX = fitInfo.numVoxels(2); nY = fitInfo.numVoxels(1); nZ = fitInfo.numVoxels(3); 
            dim = length( size(ImageOrg) );
            image2D = max( ImageOrg, [], 3); intMin = min( ImageOrg(:) ); intMax = max( ImageOrg(:) );
            
            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

            time = fitInfo.time;
            channel = fitInfo.channel;
            fitScope = fitInfo.fitScope;
            if strcmp( fitScope, 'local')
                fName = sprintf('C%d_T%d_%s_feat%d_iter%d', channel, time, fitScope, fitInfo.featureIndex, optim.iteration );
                sName = sprintf('%s_feat%d_iter%d',  fitScope, fitInfo.featureIndex, optim.iteration );
            elseif strcmp( fitScope, 'global')
                fName = sprintf('C%d_T%d_%s_iter%d', channel, time, fitScope, optim.iteration );
                sName = sprintf('%s_iter%d', fitScope, optim.iteration );
            elseif strcmp( fitScope, 'globum_add')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            elseif strcmp( fitScope, 'globum_remove')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            end
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            graphicsVerbose = 1;
            if graphicsVerbose == 1
                %  Interactive Plotting {{{
                switch state
                    case 'init'

                        % Make Figure
                        figTitle = [];
                        figName = [];
                        figPath = [];
                        f = figure( figProps{:} );
                        set( f, 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope) );
                        ax = tight_subplot(2,3, 0.05);
                        delete( ax(6) );
                        drawnow
                        pause(1)                        

                        % Original Image
                        axes( ax(1) ); %                     subplot(231)
                        img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');
                        set( ax(1), 'Tag', 'ax1');

                        % Simulated Image
    %                     subplot(232)
                        axes( ax(2) ); 
                        imageSim = FitEngine.SimulateImage( vec, fitInfo);
                        img = imagesc( max(imageSim, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated'); set( img, 'Tag', 'imgSim');
                        set( ax(2), 'Tag', 'ax2');

                        % Residual Image
    %                     subplot(233)
                        axes( ax(3) ); 
                        img = imagesc( max( abs( ImageOrg - imageSim), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Residual'); set( img, 'Tag', 'imgRes');
                        set( ax(3), 'Tag', 'ax3');

                        % Graphical Features
    %                     subplot(234)
                        axes( ax(4) ); 
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(4), 'Tag', 'ax4');

                        % Residual vs Iteration
    %                     subplot(235)
                        axes( ax(5) ); 
                        xtickformat( '%.2f'); ytickformat( '%d')
                        plotRes = plot( optim.iteration, optim.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set( plotRes,'Tag','plotRes');
                        xlabel('Iteration','interp','none'); ylabel('Resnorm','interp','none');
                        title( sprintf('Best Resnorm : %g', optim.resnorm ),'interp','none');
                        set( gca, 'FontSize', 14); grid minor; grid on
                        set( ax(5), 'Tag', 'ax5');

                        % VarX vs Iteration

                    case 'iter'
                        % Simulated Image

                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        set( f, figProps{:});
    %                     subplot(232)
                        axes( findobj('Tag', 'ax2') );
                        imageSim = FitEngine.SimulateImage( vec, fitInfo);
                        img = findobj( get( gca,'Children'), 'Tag', 'imgSim');
                        set( img, 'CData', max( imageSim, [], 3) );
                        
                        % Residual Image
    %                     subplot(233)
                        axes( findobj('Tag', 'ax3') );
                        imgRes = max( abs( ImageOrg - imageSim), [], 3); 
                        img = findobj( get( gca,'Children'), 'Tag', 'imgRes');
                        set( img, 'CData', max( imgRes, [], 3) );

                        % Graphical Features
    %                     subplot(234)
                        axes( findobj('Tag', 'ax4') );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))

                        % Update residual plot
    %                     subplot(235)
                        axes( findobj('Tag', 'ax5') );
                        plotRes = findobj( get( gca,'Children'), 'Tag', 'plotRes');
                        X = [ get( plotRes, 'Xdata') optim.iteration]; Y = [ get( plotRes, 'Ydata') optim.resnorm ];
                        set( plotRes, 'Xdata', X, 'Ydata', Y);
                        set( get( gca, 'Title'), 'String', sprintf('Best Resnorm : %g',optim.resnorm) );

                        drawnow
                        pause(0.5)

                    case 'done'
                        % No clean up tasks required for this plot function.        
                        stop = true;
                end    

                % Make a folder and save these figures
                sName = [fitInfo.saveDirectory, filesep, sName];
                export_fig( sName, '-png', '-nocrop', '-a1') 
                %  }}}

            elseif graphicsVerbose == 0

                % do all this!
                % what should i save?
                % init:
                %   time
                %   channel
                %   fitScope
                %   featureNumber
                %   featureObj
                %
                % vector labels (1st iteration)
                % vector values (all iterations) (make an array)
                % 

            end


        end
        % }}}
        % save_midfit {{{
        function stop = save_midfit( vec, optim, state, fitInfo )
            % Plots during Lsqnonlin fitting
            % Inputs :  
            %   OPTIMVALUES: properties are funccount, iteration
            %   STATE: string representing current state of optimization engine. Values are init, iter, done 
            %   STOP: A boolean to stop the algorithm.
            %   All other inputs available to error function are also available. Currently: ImageOrg( original image), fitInfo and parameters

            stop = false;

            ImageOrg = fitInfo.image;
            nX = fitInfo.numVoxels(2); nY = fitInfo.numVoxels(1); nZ = fitInfo.numVoxels(3); 
            dim = length( size(ImageOrg) );
            image2D = max( ImageOrg, [], 3); intMin = min( ImageOrg(:) ); intMax = max( ImageOrg(:) );
            
            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

            time = fitInfo.time;
            channel = fitInfo.channel;
            fitScope = fitInfo.fitScope;
            if strcmp( fitScope, 'local')
                fName = sprintf('C%d_T%d_%s_feat%d_iter%d', channel, time, fitScope, fitInfo.featureIndex, optim.iteration );
                sName = sprintf('%s_feat%d_iter%d',  fitScope, fitInfo.featureIndex, optim.iteration );
            elseif strcmp( fitScope, 'global')
                fName = sprintf('C%d_T%d_%s_iter%d', channel, time, fitScope, optim.iteration );
                sName = sprintf('%s_iter%d', fitScope, optim.iteration );
            elseif strcmp( fitScope, 'globum_add')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            elseif strcmp( fitScope, 'globum_remove')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            end
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            graphicsVerbose = 1;
            if graphicsVerbose == 1
                %  Interactive Plotting {{{
                switch state
                    case 'init'

                        % Make Figure
                        figTitle = [];
                        figName = [];
                        figPath = [];
                        f = figure( figProps{:} );
                        set( f, 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope) );
                        ax = tight_subplot(2,3, 0.05);
                        delete( ax(6) );
                        drawnow
                        pause(1)                        

                        % Original Image
                        axes( ax(1) ); %                     subplot(231)
                        img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');
                        set( ax(1), 'Tag', 'ax1');

                        % Simulated Image
    %                     subplot(232)
                        axes( ax(2) ); 
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = imagesc( max(imageSim, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated'); set( img, 'Tag', 'imgSim');
                        set( ax(2), 'Tag', 'ax2');

                        % Residual Image
    %                     subplot(233)
                        axes( ax(3) ); 
                        img = imagesc( max( abs( ImageOrg - imageSim), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Residual'); set( img, 'Tag', 'imgRes');
                        set( ax(3), 'Tag', 'ax3');

                        % Graphical Features
    %                     subplot(234)
                        axes( ax(4) ); 
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(4), 'Tag', 'ax4');

                        % Residual vs Iteration
    %                     subplot(235)
                        axes( ax(5) ); 
                        xtickformat( '%.2f'); ytickformat( '%d')
                        plotRes = plot( optim.iteration, optim.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set( plotRes,'Tag','plotRes');
                        xlabel('Iteration','interp','none'); ylabel('Resnorm','interp','none');
                        title( sprintf('Best Resnorm : %g', optim.resnorm ),'interp','none');
                        set( gca, 'FontSize', 14); grid minor; grid on
                        set( ax(5), 'Tag', 'ax5');

                        % VarX vs Iteration

                    case 'iter'
                        % Simulated Image

                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        set( f, figProps{:});
    %                     subplot(232)
                        axes( findobj('Tag', 'ax2') );
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = findobj( get( gca,'Children'), 'Tag', 'imgSim');
                        set( img, 'CData', max( imageSim, [], 3) );
                        
                        % Residual Image
    %                     subplot(233)
                        axes( findobj('Tag', 'ax3') );
                        imgRes = max( abs( ImageOrg - imageSim), [], 3); 
                        img = findobj( get( gca,'Children'), 'Tag', 'imgRes');
                        set( img, 'CData', max( imgRes, [], 3) );

                        % Graphical Features
    %                     subplot(234)
                        axes( findobj('Tag', 'ax4') );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))

                        % Update residual plot
    %                     subplot(235)
                        axes( findobj('Tag', 'ax5') );
                        plotRes = findobj( get( gca,'Children'), 'Tag', 'plotRes');
                        X = [ get( plotRes, 'Xdata') optim.iteration]; Y = [ get( plotRes, 'Ydata') optim.resnorm ];
                        set( plotRes, 'Xdata', X, 'Ydata', Y);
                        set( get( gca, 'Title'), 'String', sprintf('Best Resnorm : %g',optim.resnorm) );

                        drawnow
                        pause(0.5)

                    case 'done'
                        % No clean up tasks required for this plot function.        
                        stop = true;
                end    

                % Make a folder and save these figures
                sName = [fitInfo.saveDirectory, filesep, sName];
                export_fig( sName, '-png', '-nocrop', '-a1') 
                %  }}}

            elseif graphicsVerbose == 0

                % do all this!
                % what should i save?
                % init:
                %   time
                %   channel
                %   fitScope
                %   featureNumber
                %   featureObj
                %
                % vector labels (1st iteration)
                % vector values (all iterations) (make an array)
                % 

            end


        end
        % }}}
        % displayFinalFit {{{
        function h = displayFinalFit( Image2Fit, mainFeature, fitInfo)
            % Display Final Features and Fitting Results

            % Check if mainFeature has been dit. If not, then return out of method
            if isempty( mainFeature.featureList)
                fprintf('                       Skipping display for %s\n', mainFeature.type);  
                return
            end

            nX = size( Image2Fit, 2); nY = size( Image2Fit, 1); nZ = size( Image2Fit, 3); dim = length( size(Image2Fit) );
            image2D = max( Image2Fit, [], 3); intMin = min( Image2Fit(:) ); intMax = max( Image2Fit(:) );
            
            fName = sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope );
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';
            h = figure( figProps{:} );
            ax = tight_subplot(2, 3, 0.05);
            drawnow; pause(3)

            % Original Image
%             subplot(231)
            axes( ax(1) );
            img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');

            % Initial Simulated Image
%             subplot(232)
            axes( ax(2) );
            imageSimI = FitEngine.SimulateImage( fitInfo.fitVecs.vec, fitInfo);
            img = imagesc( max(imageSimI, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated Init'); set( img, 'Tag', 'imgSimI');

            % Final Simulated Image
%             subplot(233)
            axes( ax(3) );
            imageSimF = FitEngine.SimulateImage( fitInfo.fitResults.vfit, fitInfo);
            img = imagesc( max(imageSimF, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated Final'); set( img, 'Tag', 'imgSimF');

            % Features
%             subplot(234) 
            axes( ax(4) );
            img = imagesc( image2D ); eval( imageSets); set( img, 'Tag', 'imgFeatures'); hold on;
            mainFeature.displayFeature(gca);
            set( get(gca, 'title'), 'String', sprintf('Features : N = %d',fitInfo.featureMain.getSubFeatureNumber()))

            % Initial Residual
%             subplot(235)
            axes( ax(5) ); 
            img = imagesc( max( abs( Image2Fit - imageSimI), [], 3) ); eval( imageSets); 
            set( get(gca, 'title'), 'String', 'Image Residual Init 100%'); set( img, 'Tag', 'imgResI');

            % Final Residual
%             subplot(236)
            axes( ax(6) );
            residFrac = 100*sum( abs( Image2Fit(:) - imageSimF(:) ) ) / sum( abs( Image2Fit(:) - imageSimI(:) ) );
            img = imagesc( max( abs( Image2Fit - imageSimF), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', sprintf('Image Residual Final %.1f%%', residFrac) ); set( img, 'Tag', 'imgResF');

            sName = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ) ];
            export_fig( sName, '-png', '-nocrop', '-a1')

        end
        % }}}
        % dispImg {{{        
        function dispImg( varargin)
                
            if length(varargin) == 1 && ( isa( varargin{1}, 'double') || isa( varargin{1}, 'logical') || isa( varargin{1}, 'uint8') || isa( varargin{1}, 'single') )
                % its a single image
                img = varargin{1};
                figure; imagesc(img); axis equal; 
                if size(img,3) == 1
                    colormap gray
                end

                pos = get(gcf, 'position');
                set(gcf, 'pos', [-10000 pos(2:end)]);
                set(gcf, 'WindowState', 'maximized');
                set(gca, 'xlim', [1 size(img, 2)], 'ylim', [1 size(img, 1)], 'xtick',[], 'ytick',[] )
                
            elseif length( varargin) > 1
                % there are multiple images for comparison. The display order is
                % specified in the last argument [nrows ncols]
                
                if size( varargin{end}, 2) ~= 2
                    error('dispImg: the last argument must specify the order of plots [nRows nCols]')
                end
                
                nFig = length( varargin)-1;
                order = varargin{ end};
                nR = order(1); nC = order(2);
                img = varargin(1 : nFig);
                figure; 
                
                for jFig = 1 : nFig

                    subplot( nR, nC, jFig); imagesc( img{jFig} ); axis equal; 
                    if size( img{jFig}, 3) == 1
                        colormap gray
                    end
                    set(gca, 'xlim', [1 size(img{jFig}, 2)], 'ylim', [1 size(img{jFig}, 1)], 'xtick',[], 'ytick',[] )
                end
                pos = get(gcf, 'position');
                set(gcf, 'pos', [-10000 pos(2:end)]);
                set(gcf, 'WindowState', 'maximized');
                    
            else
                error('Oops, something went wrong!')
            end

        end
        % }}}

        % }}}
    
        % saveFinalFit {{{
        function h = saveFinalFit( Image2Fit, mainFeature, fitInfo)
            % Save data for this particular fit

            % Book-keeping
            channel = uint8(fitInfo.channel);
            time = uint16(fitInfo.time);
            fitScope = fitInfo.fitScope;
%             Image2Fit = uint16 (Image2Fit);

            % Images Simulated
            if ~isempty( mainFeature.featureList) && Cell.checkFrameViability( Image2Fit, time)
                imageSimI = uint16( FitEngine.SimulateImage( fitInfo.fitVecs.vec, fitInfo) );
                imageSimF = uint16( FitEngine.SimulateImage( fitInfo.fitResults.vfit, fitInfo) );
            else
                imageSimI = [];
                imageSimF = [];
            end

            % Main feature
            featureMainStruct = fitInfo.featureMain.saveAsStruct();

            f_data = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ), '.mat' ];
            save( f_data, 'channel', 'time', 'fitScope', 'imageSimI', 'imageSimF', 'featureMainStruct', '-v7');

        end
        % }}}
        
        % PrintInfo {{{
        function PrintInfo( obj)
            % Prints information about the cell

            fprintf('-----------------------   C E L L    I N F O   ------------------------\n');
            fprintf('-----------------------------------------------------------------------\n\n');
            fprintf( 'Type : %s\n', obj.type)
            fprintf( 'ImageData : \n')
            ImageData.PrintInfo( obj.imageData);
            disp( ['Features = ', strjoin( obj.featuresInChannels, ' - ') ] ) 
            fprintf('----------------------------------------------------------------------\n');

        end
        % }}}

    end

    methods ( Abstract = true )
    
        obj = EstimateFeatures( obj, image, time, channelTrue, channelIdx)

    end

end
