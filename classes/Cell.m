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
        featureProps = Cell.GetFeatureProps();
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
                
                if ~isfield( obj.params, 'newEstimateEveryT')
                    obj.params.newEstimateEveryT = 0;
                end
                
                if obj.params.timeReversal && obj.params.newEstimateEveryT
                    obj.params.timeReversal = 0;
                end
                if obj.params.timeReversal
                    obj.params.timeReversal = 0;
                    obj.FindFeaturesChannel( jChannel);
                    obj.params.timeReversal = 1;
                else
                    obj.FindFeaturesChannel( jChannel);
                end
                
                % Reverse time
                if obj.params.timeReversal
                    fprintf('Time reversed!\n')
                    % Basic cleanup
                    obj.featureList = cell( length( obj.channelsToFit), size( obj.imageData.GetImage, 4) );
                    obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
                    % fit again
                    obj.FindFeaturesChannel( jChannel);
                end

            end

        end
        % }}}

        % FindFeaturesChannel {{{
        function obj = FindFeaturesChannel( obj, jChannel)

            % get actual channel number based on channel index
            cChannel = obj.channelsToFit( jChannel);

            % Get lifetime for this cell
            lifetime = obj.imageData.GetLifetime;
            lifetimes = lifetime(1):lifetime(2);
            
            % Reverse time order if parameter specified
            if obj.params.timeReversal
                lifetimes = flip( lifetimes);
            end

            % Check for existence of images in all the frames and duplicate images if necessary
            [obj, imgNew] = CheckUpdateImages( obj, cChannel);

            for jTime = 1: length(lifetimes)
                
                cTime = lifetimes(jTime);

                fprintf( 'C-%d Time = %d\n', jChannel, cTime) 

                % Initialize for this CT step
                params = obj.params;
                params.channelTrue = cChannel;
                params.channelIdx = jChannel;
                params.time = cTime;
                params.timeIdx = jTime;
                params.saveDirectory = [ obj.params.saveDirectory , filesep, sprintf( 'C%d_T%d', cChannel, cTime) ];
                mkdir( params.saveDirectory);

                % Find features for this CT frame
                obj.FindFeaturesFrame( imgNew, params) 

            end
            disp('-----------------------------------------------------------------------')
            
        end
        % }}}
        
        % FindFeaturesFrame {{{
        function FindFeaturesFrame( obj, Image, parameters)
            % Find and fit features for a single CT frame
            
            plot_estimate_and_skip = 1;
            
            % Decide if fit should be performed. Skip otherwise
            [obj, status] = DecideToFit( obj, Image, parameters);
            if ~status
                return
            end

            % Get the image for this frame
            Image2Fit = Image(:,:,:, parameters.time, parameters.channelTrue);
            
            % Estimate the features based on an estimation routine (defined in specialized sub-class )
            disp('Estimating features...')
            % Good estimation is key to good optimization in low SnR images
            obj.EstimateFeatures( Image2Fit, parameters.time, parameters.channelTrue, parameters.channelIdx,obj.params.timeReversal, obj.params.newEstimateEveryT);
            mainFeature = obj.featureList{ parameters.channelIdx , parameters.time};
            obj.syncFeatureMap( parameters.channelIdx, parameters.time);
            
            if plot_estimate_and_skip
                nX = size( Image2Fit, 2); nY = size( Image2Fit, 1); nZ = size( Image2Fit, 3); dim = length( size(Image2Fit) );
                image2D = max( Image2Fit, [], 3); intMin = min( Image2Fit(:) ); intMax = max( Image2Fit(:) );
                imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

                figure;
                ax = tight_subplot(1, 2, 0.05);
                axes( ax(1) );
                img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original');
                axes( ax(2) );
                img = imagesc( image2D ); eval( imageSets); hold on;
                mainFeature.displayFeature(gca); set( get(gca, 'title'), 'String', 'Features');
                suptitle( num2str(parameters.time))
                drawnow; pause(1);close all;
                return
            end
            
            % Pre-fit adjustments
            mainFeature.preOptimize();
            
            % Prepare fit params
            params = obj.params.fit;
            params.channel = parameters.channelTrue;
            params.time = parameters.time;
            params.saveDirectory = parameters.saveDirectory;
            params.timeReversal = obj.params.timeReversal;
            
            % Fit Features via Fit Engine.
            fitEngine = FitEngine( Image2Fit, mainFeature, params);
            fitEngine = fitEngine.Optimize();
            obj.featureList{ parameters.channelIdx , parameters.time} = fitEngine.GetFeature();
            obj.featureList{ parameters.channelIdx , parameters.time}.finalizeAddedFeatures();
            
            obj.syncFeatureMap( parameters.channelIdx, parameters.time);
            close all

        end
        % }}}
        
        % CheckUpdateImages {{{
        function [obj, imgNew] = CheckUpdateImages( obj, channel)

            % Get Image
            imgOld = obj.imageData.GetImage();
            imgNew = 0*imgOld;

            % Get lifetime
            lifetime = obj.imageData.GetLifetime();
            life_vec = lifetime(1) : lifetime(2);

            % Find good frames
            frames_good = Cell.FindGoodFrames( imgOld(:,:,:,:,channel), life_vec);
    
            % Print good frames
%             disp( 'Printing good frames...')
%             disp( life_vec( find(frames_good) ) )

            if sum( frames_good) == 0
                error('no good frames over here. Something must be wrong...')
            end

            % Loop over time and select the new image for bad frames using good frames
            for jt = 1: length( life_vec )

                % if frame is good, do nothing
                if frames_good( jt)
                    imgNew(:,:,:, life_vec( jt), channel) = imgOld(:,:,:, life_vec( jt), channel);
                    continue

                % If frame is bad, then update the image to one of a good frame
                else
                    % find next good frame
                    jt_next = (jt-1) + find( frames_good( jt: end), 1, 'first');

                    % find prev good frame
                    jt_prev = find( frames_good( 1: jt), 1, 'last');

                    % If no prev good frame, then set image to the next good frame
                    if isempty( jt_prev)
                        imgNew(:,:,:, life_vec( jt), channel) = imgOld(:,:,:, life_vec( jt_next), channel);
                        continue
                    % If no next good frame, then set image to the prev good frame
                    elseif isempty( jt_next) 
                        imgNew(:,:,:, life_vec( jt), channel) = imgOld(:,:,:, life_vec( jt_prev), channel);
                        continue
                    end

                    % If both prev and next good frames exist, use the one closest in time 
                    if abs(jt-jt_prev) < abs(jt-jt_next)
                        imgNew(:,:,:, life_vec( jt), channel) = imgOld(:,:,:, life_vec( jt_prev), channel);
                    else
                        imgNew(:,:,:, life_vec( jt), channel) = imgOld(:,:,:, life_vec( jt_next), channel);
                    end
                end
            end
        
        end
        % }}}

        % DecideToFit {{{
        function [obj, status] = DecideToFit( obj, Image, parameters)
            % In the case of missing frames, fitting must be skipped and features must be copied.

            status = 1;
            channel = parameters.channelTrue;
            channelIdx = parameters.channelIdx;
            time = parameters.time;
            
            % Get lifetimes
            lifetime = obj.imageData.GetLifetime();
            life_vec = lifetime(1): lifetime(2);

            % Check if this is the first frame:
            % If yes, then fit it
            % If not first frame:
            %   Then check if frame is the same as last frame. If it is, then do the following:
            %       1. Skip fitting
            %       2. Copy features from old frame to current frame
            %       3. Save the final fit
            
            currFrame = Image(:,:,:,time,channel);
            if time == 1
                prevFrame = 0*currFrame;
            else
                prevFrame = Image(:,:,:,time-1,channel);
            end

            if time ~= life_vec(1) && all( currFrame(:) == prevFrame(:) )

                % Skip fitting
                status = 0;

                % Copy features from the previous frame
                obj.featureList{ channelIdx, time} = obj.featureList{ channelIdx, time-1}.copyDeep();

                % Prepare to save the fit
                p.channel = channel;
                p.time = time; 
                p.fitScope = 'final'; 
                p.saveDirectory = parameters.saveDirectory;
                p.featureMain = obj.featureList{ channelIdx, time};

                % Save the final fit
                Cell.saveFinalFit( Image, obj.featureList{ channelIdx, time}, p);
                    
            end
        end
        % }}}
        
        % PropagateOldFeature {{{
        function obj = PropagateOldFeature(obj, idxChannel, cChannel, cTime, timeReverse)

            disp('- Propagating old feature') 

            % Find the most recent good frame
            if timeReverse
                bestFrame = cTime+1;
            else
                bestFrame = cTime-1;
            end
            while obj.featureList{ idxChannel, bestFrame} ~= obj.featureList{ idxChannel, bestFrame}
                if timeReverse
                    bestFrame = cTime+1;
                else
                    bestFrame = cTime-1;
                end
            end

            % Duplicate feature from best recent frame
            obj.featureList{ idxChannel, cTime} = obj.featureList{ idxChannel, bestFrame}.copyDeep();

            % Ensure the image is from the actual frame
            Image = obj.imageData.GetImage();
            obj.featureList{ idxChannel, cTime}.image = im2double( Image(:,:,:, cTime, cChannel) );

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
                amplitude = interp2( imageIn, xq, yq, 'linear');

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
        % findAmplitudeAlongCurve {{{
        function amplitude = findAmplitudeAlongCurve( imageIn, coeffs)
            % used to measure intensity while using polynomial coefficients
            
            dim = length (size( imageIn) );
            t = linspace(0,1);
            px = coeffs{1};
            py = coeffs{2};
            x = polyval(px, t);
            y = polyval(py, t);
            if dim == 3
                pz = coeffs{3};
                z = polyval(pz, t);
            end
            if dim ==2
                coord = round( [x; y]);
            elseif dim ==3
                coord = round( [x; y; z]);
            end
            % keep unique coordinates
            idxRm = [];
            for jCrd = 1 : size(coord, 2)-1
                if dim ==2
                    if ( coord(1, jCrd) == coord(1, jCrd+1) ) && ( coord(2, jCrd) == coord(2, jCrd+1) )
                        idxRm = [idxRm, jCrd];
                    end
                elseif dim ==3
                    if ( coord(1, jCrd) == coord(1, jCrd+1) ) && ( coord(2, jCrd) == coord(2, jCrd+1) ) && ( coord(3, jCrd) == coord(3, jCrd+1) )
                        idxRm = [idxRm, jCrd];
                    end
                end
            end

            coord(:, idxRm) = [];
            if dim ==2
                idx = sub2ind( size(imageIn), coord(2, :), coord(1, :) );
            elseif dim ==3
                idx = sub2ind( size(imageIn), coord(2, :), coord(1, :), coord(3, :) );
            end
            amplitude = imageIn( idx);
            
        end
        % }}}
        % findAmplitudeAlongCurveCoords {{{
        function amplitude = findAmplitudeAlongCurveCoords( imageIn, coord)
            % used to measure intensity while using coordinates
            
            dim = length (size( imageIn) );
            
            if dim ==2
                idx = sub2ind( size(imageIn), coord(2, :), coord(1, :) );
            elseif dim ==3
                idx = sub2ind( size(imageIn), coord(2, :), coord(1, :), coord(3, :) );
            end
            amplitude = imageIn( idx);
            
        end
        % }}}
        % radIntegrate2D {{{
        function [phiIntensity, phiValues] = radIntegrate2D( imageIn, startPoint, rmin, rmax)

            if length( size(imageIn) ) ~= 2
                error( 'radIntegrate2D : input image must be 2-dimensional')
            end
            
            if nargin < 4
                rmax = 15;
            end
            if nargin < 3
                rmin = 1;
            end

            % Radially integrate around startPoint 
            PhiStep = deg2rad(2);
            numVoxelsX = size( imageIn, 2);
            numVoxelsY = size( imageIn, 1);
            x0 = startPoint(1);
            y0 = startPoint(2);

            % Define the step sizes to use for the angular sweep in phi. 
            PhiStep = deg2rad(2);
%             rmin = 5; % 1
%             rmax = 15; % 15
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
        function [radIntensityMean, radValues] = phiIntegrate2D( imageIn, startPoint, rmax)

            if length( size(imageIn) ) ~= 2
                error( 'radIntegrate2D : input image must be 2-dimensional')
            end

            ignoreZerosInMean = 1;

            % Radially integrate around startPoint 
            radStep = 0.25;
            rmin = 1;
            if nargin < 3
                rmax = 100;
            end
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
            Xcirc = x0 + radValues'.*cos( phiValues);
            Ycirc = y0 + radValues'.*sin( phiValues);

            % Get intensity values of these scattered points in 2D space
            Intcirc = interp2( im2double( imageIn), Xcirc, Ycirc, 'linear', 0);
           
            % Sum over phi to get a vector of Intensity values at varying radia
            radIntensity = sum( Intcirc, 2)';

            % Find a mean
            if ignoreZerosInMean
                radIntensityMean = radIntensity ./ sum( logical( Intcirc), 2)';
            else
                radIntensityMean = radIntensity / size( Intcirc,2); 
            end

        end
        % }}}
        % checkFrameViability {{{
        function status = checkFrameViability( Image2Fit, jTime)

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
        % FindGoodFrames {{{
        function frames_good = FindGoodFrames( Image, lifetimes)

            % find all good frames
            frames_good = 0*lifetimes;

            for j = 1: length(lifetimes)
                
                jt = lifetimes( j);
                Img = Image(:,:,:, jt);
                if max( Img(:)) ~= min( Img(:)) 
                    frames_good( j) = 1;
                end
            end
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

            Lmin = 4.5;

            % This will look at radial amplitude alon the Line to find when the intensity falls below a threshold
            threshIntensity = median( image2D( image2D ~= 0) ) + 1*std( image2D( image2D ~= 0) );
%             median( image2D( image2D ~= 0) ) 
%             std( image2D( image2D ~= 0) )
%             error('stop here')

            % Define what radia to use:
            rmin = 1;
            rmax = 40;
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
            radialIntensity = imgaussfilt( radialIntensity, 3);

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
            nX = fitInfo.numVoxels(2); nY = fitInfo.numVoxels(1); 
            try
                nZ = fitInfo.numVoxels(3);
            end
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
            if fitInfo.timeReversal
                sName = [fitInfo.saveDirectory, filesep, sName, '_reverse'];
            else
                sName = [fitInfo.saveDirectory, filesep, sName];
            end
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 2000 1280 720]};

            graphicsVerbose = 1;
            if graphicsVerbose == 1
                %  Interactive Plotting {{{
                switch state
                    case 'init'

                        % Make Figure
                        figTitle = []; figName = []; figPath = [];
                        f = figure( figProps{:} );
                        set( f, 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope) );
                        ax = tight_subplot(2,3, 0.05); drawnow; pause(0.2)

                        % Original Image
                        axes( ax(1) );
                        img = imshow( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');
                        set( ax(1), 'Tag', 'ax1');

                        % Simulated Image
                        axes( ax(2) ); 
                        imageSim = FitEngine.SimulateImage( vec, fitInfo);
                        img = imagesc( max(imageSim, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated'); set( img, 'Tag', 'imgSim');
                        set( ax(2), 'Tag', 'ax2');

                        % Residual Image
                        axes( ax(3) ); 
                        img = imagesc( max( abs( ImageOrg - imageSim), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Residual'); set( img, 'Tag', 'imgRes');
                        set( ax(3), 'Tag', 'ax3');

                        % Graphical Features
                        axes( ax(4) );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(4), 'Tag', 'ax4');

                        % Vector features Z
                        axes( ax(5) ); hold on; ax(5).Color = 'k';
                        axis equal; axis ij; set( gca, 'xlim', [1 nX], 'ylim', [1 nY], 'XTick', [], 'YTick', [], 'FontSize', 14)
                        colormap( ax(5), cool);
                        % boundary
                        [B,~] = bwboundaries(image2D > 0,'noholes');
                        plot( B{1}(:,2), B{1}(:,1), 'color', 'w', 'linewidth',3)
                        if strcmp( fitScope, 'local')
                            fitInfo.featureCurrent.displayFeature(gca, nZ);
                        else
                            fitInfo.featureCurrent.displayFeature(gca, 1);
                        end
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(5), 'Tag', 'ax5'); 

                        % Residual vs Iteration
                        axes( ax(6) ); 
                        xtickformat( '%.2f'); ytickformat( '%d')
                        plotRes = plot( optim.iteration, optim.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set( plotRes,'Tag','plotRes');
                        xlabel('Iteration','interp','none'); ylabel('Resnorm','interp','none');
                        title( sprintf('Best Resnorm : %g', optim.resnorm ),'interp','none');
                        set( gca, 'FontSize', 14); grid minor; grid on
                        set( ax(6), 'Tag', 'ax6');

                    case 'iter'
                        % Simulated Image
                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        set( f, figProps{:});
                        axes( findobj('Tag', 'ax2') );
                        imageSim = FitEngine.SimulateImage( vec, fitInfo);
                        img = findobj( get( gca,'Children'), 'Tag', 'imgSim');
                        set( img, 'CData', max( imageSim, [], 3) );
                        
                        % Residual Image
                        axes( findobj('Tag', 'ax3') );
                        imgRes = max( abs( ImageOrg - imageSim), [], 3); 
                        img = findobj( get( gca,'Children'), 'Tag', 'imgRes');
                        set( img, 'CData', max( imgRes, [], 3) );

                        % Graphical Features
                        axes( findobj('Tag', 'ax4') );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))

                        % Update residual plot
                        axes( findobj('Tag', 'ax6') );
                        plotRes = findobj( get( gca,'Children'), 'Tag', 'plotRes');
                        X = [ get( plotRes, 'Xdata') optim.iteration]; Y = [ get( plotRes, 'Ydata') optim.resnorm ];
                        set( plotRes, 'Xdata', X, 'Ydata', Y);
                        set( get( gca, 'Title'), 'String', sprintf('Best Resnorm : %g',optim.resnorm) );

                        % Vector features Z
                        axes( findobj('Tag', 'ax5') ); hold on;
                        cla( findobj('Tag', 'ax5')); 
                        colormap( findobj('Tag', 'ax5'), cool);
                        % boundary
                        [B,~] = bwboundaries(image2D > 0,'noholes');
                        plot( B{1}(:,2), B{1}(:,1), 'color', 'w', 'linewidth',3)
                        if strcmp( fitScope, 'local')
                            fitInfo.featureCurrent.displayFeature(gca, nZ);
                        else
                            fitInfo.featureCurrent.displayFeature(gca, 1);
                        end
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        drawnow; pause(0.2)

                    case 'done'
                        stop = true;
                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        close(f)
                end    

                % Make a folder and save these figures
                try
                    export_fig( sName, '-png', '-nocrop', '-a1') 
                catch
                    disp('fig not saved in Cell.plot_midfit')
                end
                %  }}}
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
            imageSimI = FitEngine.SimulateImage( fitInfo.fitVecOld, fitInfo.fitInfoOld);
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
        % saveFinalFit {{{
        function h = saveFinalFit( Image2Fit, mainFeature, fitInfo)
            % Save data for this particular fit
            
            % Book-keeping
            channel = uint8(fitInfo.channel);
            time = uint16(fitInfo.time);
            fitScope = fitInfo.fitScope;
%             Image2Fit = uint16 (Image2Fit);

            % Images Simulated
            try
                imageSimI = fitInfo.imageSimI;
                imageSimF = fitInfo.imageSimF;
            catch
                imageSimI = [];
                imageSimF = [];
            end
            
            % Main feature
            featureMainStruct = fitInfo.featureMain.saveAsStruct();

            if fitInfo.timeReversal
                f_data = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ), '_reverse.mat' ];
            else
                f_data = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ), '.mat' ];
            end
            save( f_data, 'channel', 'time', 'fitScope', 'imageSimI', 'imageSimF', 'featureMainStruct', '-v7');

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

        % GetFeatureProps {{{
        function props = GetFeatureProps()
        % This is a structure of fit and graphics properties for all features.

            % Environment {{{
            env1 = {'background'};
            env2 = {'background', 'backgroundNuclear'};
            % }}}
            
            % Basic Elements {{{
            % Spot
            spot.fit{2} = {'position', 'amplitude', 'sigma'};
            spot.fit{3} = {'position', 'amplitude', 'sigma'};
            spot.graphics.magenta = {'Color', [0.7 0 0.7] , 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2};
            spot.graphics.green = {'Color', [0 1 0] , 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2};
            spot.graphics.blue = {'Color', [0 0 1] , 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2};
            spot.graphics.red = {'Color', [1 0 0] , 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2};
            spot.graphics.purple = {'Color', [0.35 0.06 0.51] , 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 5};

            % Line
%             line.fit{2} = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
%             line.fit{3} = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            line.fit{2} = {'startPosition', 'length','theta', 'amplitude', 'sigma'};
            line.fit{3} = {'startPosition', 'length','theta', 'amplitude', 'sigma'};
            line.graphics.magenta = {'Color', [0.7 0 0.7] , 'LineWidth', 2};
            line.graphics.green = {'Color', [0 1 0] , 'LineWidth', 2};
            line.graphics.blue = {'Color', [0 0 1] , 'LineWidth', 2};
            line.graphics.red = {'Color', [1 0 0] , 'LineWidth', 2};
            line.graphics.redWide = {'Color', [1 0 0] , 'LineWidth', 4};
            line.graphics.greenWide = {'Color', [0 1 0] , 'LineWidth', 4};
            
            % Curve
            curve.fit{2} = {'startPosition','cX','cY','amplitude','sigma'};
            curve.fit{3} = {'startPosition', 'cX','cY','cZ', 'amplitude', 'sigma'};
            curveMT.fit{2} = {'origin', 'thetaInit', 'normalVec', 'amplitude', 'L','sigma'};
            curveMT.fit{3} = {'origin', 'thetaInit', 'normalVec', 'amplitude', 'L','sigma'};
            curve.graphics.magenta = {'Color', [0.7 0 0.7] , 'LineWidth', 2};
            curve.graphics.green = {'Color', [0 1 0] , 'LineWidth', 2};
            curve.graphics.blue = {'Color', [0 0 1] , 'LineWidth', 2};
            curve.graphics.red = {'Color', [1 0 0] , 'LineWidth', 2};

            % Bundle 
%             bundle.fit{2} = {'cX','cY', 'amplitude', 'sigma', 'T', 'ef'};
%             bundle.fit{3} = {'cX','cY','cZ', 'amplitude', 'sigma', 'T', 'ef'};
            bundle.fit{2} = {'origin', 'thetaInit', 'normalVec', 'amplitude', 'T', 'L','ef', 'sigma'};
            bundle.fit{3} = {'origin', 'thetaInit', 'normalVec', 'amplitude', 'T', 'L','ef', 'sigma'};
            bundle.graphics.magenta = {'Color', [0.7 0 0.7] , 'LineWidth', 2};
            bundle.graphics.green = {'Color', [0 1 0] , 'LineWidth', 2};
            bundle.graphics.blue = {'Color', [0 0 1] , 'LineWidth', 2};
            bundle.graphics.red = {'Color', [1 0 0] , 'LineWidth', 2};

            props.spot = spot;
            props.line = line;
            props.curve = curve;
            props.bundle = bundle;
            % }}}
            
            % Organizers {{{
            % Line Aster
            asterLine.fit{2}.spot = spot.fit{2};
            asterLine.fit{3}.spot = spot.fit{3};
            asterLine.fit{2}.line = line.fit{2}(2:end);
            asterLine.fit{3}.line = line.fit{3}(2:end);
            asterLine.graphics.spot = spot.graphics.blue;
            asterLine.graphics.line = line.graphics.green;
            props.asterLine = asterLine;

            % Curve Aster
            asterCurve.fit{2}.curve = curveMT.fit{2};
            asterCurve.fit{3}.curve = curveMT.fit{3};
            asterCurve.graphics.spot = spot.graphics.blue;
            asterCurve.graphics.curve = bundle.graphics.red;
            props.asterCurve = asterCurve;
            % }}}

            % Master Organizers {{{
            % Spindle
            spindle.fit{2}.line = line.fit{2};
            spindle.fit{2}.curve = curveMT.fit{2};
            spindle.fit{2}.aster = asterLine.fit{2};
            spindle.fit{2}.aster.spot(1) = [];
            spindle.fit{2}.Environment = env1;
            spindle.fit{3}.line = line.fit{3};
            spindle.fit{3}.curve = curveMT.fit{3};
            spindle.fit{3}.aster = asterLine.fit{3};
            spindle.fit{3}.aster.spot(1) = [];
            spindle.fit{3}.Environment = env1;
            spindle.graphics.line = line.graphics.greenWide;
            spindle.graphics.curve = curve.graphics.red;
            spindle.graphics.aster.spot = spot.graphics.purple;
            spindle.graphics.aster.line = line.graphics.green;
            spindle.fit{2}.spot = spot.fit{2};
            spindle.fit{3}.spot = spot.fit{3};
            spindle.graphics.spot = spot.graphics.purple;
            props.spindle = spindle;

            % Monopolar Aster
            monopolarAster.fit{2}.aster = asterLine.fit{2};
            monopolarAster.fit{3}.aster = asterLine.fit{3};
            monopolarAster.fit{2}.Environment = env1;
            monopolarAster.fit{3}.Environment = env1;
            monopolarAster.graphics.aster.spot = spot.graphics.purple;
            monopolarAster.graphics.aster.line = line.graphics.green;
            monopolarAster.fit{2}.spot = spot.fit{2};
            monopolarAster.fit{3}.spot = spot.fit{3};
            monopolarAster.graphics.spot = curve.graphics.blue;
            monopolarAster.fit{2}.line = line.fit{2};
            monopolarAster.fit{3}.line = line.fit{3};
            monopolarAster.graphics.line = curve.graphics.red;
            props.monopolarAster = monopolarAster;

            % Interphase Bank
            intBank.fit{2}.aster = asterCurve.fit{2};
            intBank.fit{3}.aster = asterCurve.fit{3};
            intBank.fit{2}.Environment = env1;
            intBank.fit{3}.Environment = env1;
            intBank.graphics.aster.spot = spot.graphics.blue;
            intBank.graphics.aster.curve = curve.graphics.red;
            intBank.fit{2}.spot = spot.fit{2};
            intBank.fit{3}.spot = spot.fit{3};
            intBank.graphics.spot = spot.graphics.blue;
            intBank.fit{2}.curve = bundle.fit{2};
            intBank.fit{3}.curve = bundle.fit{3};
            intBank.graphics.curve = curve.graphics.red;
            props.intBank = intBank;
            % }}}

        end
        
        function routine = SetFitRoutines()
            
            % Routine type:
            %   1) L : local
            %   2) G : global
            %   3) LG: local+global
            %   4) GA: global+Add features
            %   5) def: use default parameters
            % routine 1
            % background
            routine{1}.type = 'G';
            routine{1}.props = {};
            
            % routine 2
            % origin, thetaInit, normalVec, L
            routine{2}.type = 'L';
            routine{2}.props = {'origin', 'thetaInit', 'normalVec','L'};
            
            % routine 3
            % origin, thetaInit, normalVec, L
            routine{3}.type = 'L';
            routine{3}.props = {'amplitude', 'ef', 'T', 'sigma'};
            
            % routine 4
            % origin, thetaInit, normalVec, L
            routine{4}.type = 'GA';
            routine{4}.props = {'normalVec', 'L','amplitude', 'T','ef'};
            
        end

    end
    % }}}

    methods ( Abstract = true )
    
        obj = EstimateFeatures( obj, image, time, channelTrue, channelIdx)

    end

end
