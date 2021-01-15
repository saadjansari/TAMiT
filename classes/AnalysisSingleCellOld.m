classdef AnalysisSingleCellOld < handle
    % This is the data analysis class. It is basically a container for all activities analysis related
    properties
        path
        cellType
        features
        channels = []
        times = []
        timeStep
        sizeVoxels = [0.1067 0.1067 0.5]; %default
        numVoxels 
        data
        folderName
        goodFrame = 1
        flagMovie
        flagGraph
        simImageMT
        fwdreverse = 0
    end

    methods (Access = public )

        % AnalysisSingleCell {{{
        function obj = AnalysisSingleCell( pathResCell, cellType, channels, features, timeStep, sizeVoxels, flagMovie, flagGraph)
            % Contruct the Analysis Object
            
            obj.path = pathResCell;
            if ~strcmp(cellType, 'Mitosis') && ~strcmp( cellType, 'Interphase') && ~strcmp( cellType, 'Monopolar') && ~strcmp( cellType, 'MitosisBud')
                error('AnalysisBank: unknown cellType, allowed values are Mitosis and Interphase and Monopolar.')
            end

            obj.cellType = cellType;

            obj.channels = channels;

            if ~any( strcmp(features, 'Microtubule') ) && ~any( strcmp( features, 'Kinetochore') ) && ~any( strcmp( feature, 'Cut7' ) )
                error('AnalysisBank: unknown feature, allowed values are Microtubule, Kinetochore, and Cut7.')
            end

            obj.features = features;

            if nargin > 4
                obj.timeStep = timeStep;
            end

            if nargin > 5
                obj.sizeVoxels = sizeVoxels;
            end

            [~, obj.folderName, ~] = fileparts( pathResCell);

            if nargin > 6
                obj.flagMovie = flagMovie;
                obj.flagGraph = flagGraph;
            end

        end
        % }}}
        
        % Analyze {{{
        function Analyze( obj)
            % Analyze the object

            % Analyze each channel
            for jChannel = 1: length( obj.channels)

                % Find times for analysis
                obj.times = obj.findTimes( obj.channels( jChannel));
            
                % Analysis
                fname = [obj.path,filesep,'analysisData_',num2str(jChannel),'.mat'];
                disp( 'temp bypass of analysis file check')
                fprintf('   Analysis: channel %d with feature %s...\n', obj.channels( jChannel), obj.features{ jChannel}); 
                dat = obj.analyzeChannel( jChannel); save( fname, 'dat', '-v7.3');
%                 if exist(fname) ~= 2
%                     fprintf('   Analysis: channel %d with feature %s...\n', obj.channels( jChannel), obj.features{ jChannel}); 
%                     dat = obj.analyzeChannel( jChannel); save( fname, 'dat', '-v7.3'); 
%                 else
%                     fprintf('   Analysis: data exists, loading from file...\n')
%                     load( fname); obj.data{jChannel} = dat;
%                 end

                % graph channel analysis
                if obj.flagGraph
                    fprintf('   Graphing...\n')
                    obj.graphChannel( jChannel);
                end

                % Make movie
                if obj.flagMovie 
                    fprintf('   Directing movies...\n')
                    % obj.temp_fig1( jChannel);
                    % obj.temp_fig3( jChannel);
%                     obj.makeMovie( jChannel);
                    if obj.fwdreverse
                        obj.makeMovieFwdReverse( jChannel);
                    end
                end
                
                % Tracking
                obj.trackChannel( jChannel);

            end

        end
        % }}}

        % analyzeChannel{{{
        function data = analyzeChannel( obj, jChannel)
            % Analyze this frame
            
            global COUNTER
            COUNTER =1;
            
            for jFrame = 1 : length( obj.times)
                fprintf('      Analyzing time %d...\n', obj.times( jFrame)); 
                obj.analyzeFrame( jChannel, jFrame);
            end
            obj.data{jChannel}.tag = obj.features{jChannel};
            data = obj.data{jChannel};
            
        end
        % }}}

        % analyzeFrame {{{
        function analyzeFrame( obj, jChannel, jTime)
            % Analyze the time given in cTime

            global COUNTER
            
            % Current channel 
            channel = obj.channels( jChannel);
            time = obj.times( jTime);

            % Load the data file if available
            subfold = ['C' num2str( channel) '_T' num2str( time)]; 
            timeDataFile = [obj.path, filesep, subfold, filesep, subfold, '_final.mat']; 
            timeData = load( timeDataFile);
            obj.goodFrame = jTime;

            % Initialize main feature from loaded data
            eval( ['mainFeatureFwd = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );
            
            % Try loading reverse time data
            try
%                 disp('loading reverse file')
                % Load the data file if available
                timeDataFile = [obj.path, filesep, subfold, filesep, subfold, '_final_reverse.mat']; 
                timeData = load( timeDataFile);
                obj.fwdreverse = 1;
                eval( ['mainFeatureRev = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );

                mainFeatureRev.ID = COUNTER;
                mainFeatureRev.updateFeatureIDs();
                mainFeatureRev.updateFeatureMap();
                obj.data{jChannel}.featuresRev{jTime} = mainFeatureRev;
                mainFeatureFwd.ID = COUNTER;
                mainFeatureFwd.updateFeatureIDs();
                mainFeatureFwd.updateFeatureMap();
                obj.data{jChannel}.featuresFwd{jTime} = mainFeatureFwd;
                
                % Forward-Reverse linking
                mainFeature = obj.linkForwardReverse( mainFeatureFwd, mainFeatureRev);
            catch
                mainFeatureRev = mainFeatureFwd.copyDeep();
                mainFeatureRev.ID = COUNTER;
                mainFeatureRev.updateFeatureIDs();
                mainFeatureRev.updateFeatureMap();
                obj.data{jChannel}.featuresRev{jTime} = mainFeatureRev;
                mainFeatureFwd.ID = COUNTER;
                mainFeatureFwd.updateFeatureIDs();
                mainFeatureFwd.updateFeatureMap();
                obj.data{jChannel}.featuresFwd{jTime} = mainFeatureFwd;
                mainFeature = mainFeatureFwd.copyDeep();
            end
            mainFeature.ID = COUNTER;
            mainFeature.updateFeatureIDs();
            mainFeature.updateFeatureMap();
            
            % Run Feature Specific analysis
            obj.analyzeFeature( mainFeature, jChannel, jTime);
            
            switch mainFeature.type
                case 'SpindleNew'
                    for ja = 2 : 3
                        for jc = 2 : length(mainFeature.featureList{ja}.featureList)
                            mainFeature.featureList{ja}.featureList{jc}.display = {'Color', [1 0 0] , 'LineWidth', 2};
                        end
                    end
                case 'Spindle'
                    for ja = 2 : 3
                        for jc = 2 : length(mainFeature.featureList{ja}.featureList)
                            mainFeature.featureList{ja}.featureList{jc}.display = {'Color', [0 1 0] , 'LineWidth', 2};
                        end
                    end
                case 'MonopolarAster'
                    for jc = 2 : length(mainFeature.featureList{1}.featureList)
                        mainFeature.featureList{1}.featureList{jc}.display = {'Color', [0 1 0] , 'LineWidth', 2};
                    end
            end

            % Save information for making movies
            obj.data{jChannel}.features{jTime} = mainFeature;

        end
        % }}}

        % analyzeFeature {{{
        function analyzeFeature( obj, mainFeature, jChannel, jTime)
            % Run analysis on the main feature by deciding what type of feature exists

            switch obj.features{jChannel}
                case 'Microtubule'
                    obj.analyzeFeatureMicrotubule( mainFeature, jChannel, jTime);
                case 'Kinetochore'
                    obj.analyzeFeatureKinetochore( mainFeature, jChannel, jTime);
                case 'Cut7'
                    obj.analyzeFeatureCut7( mainFeature, jChannel, jTime);
                otherwise
                    error('analyzeFrame: unknown feature')
            end

        end
        % }}}
        
        % Microtubules {{{
        % analyzeFeatureMicrotubule {{{
        function analyzeFeatureMicrotubule( obj, mainFeature, jChannel, jTime)
            % Analyze some type of microtubule feature (e.g. Spindle, Monopolar, InterphaseBank)
            
            % Simulate the mt image
            mainFeature.fillParams( size(mainFeature.image) );
            obj.simImageMT(:,:,:,jTime) = mainFeature.simulateFeature( size(mainFeature.image) );

            switch mainFeature.type
                case 'Spindle'

                    % Spindle Length
                    %[ obj.data{jChannel}.spindleLength(jTime) ] = obj.analyzeSpindle( mainFeature, jTime);
                
                case 'SpindleNew'

                    % Spindle Length
                    %[ obj.data{jChannel}.spindleLength(jTime) ] = obj.analyzeSpindle( mainFeature, jTime);

                case 'MonopolarAster'

                    % Pole position, number of microtubules and tip distance
                    %[ obj.data{jChannel}.pole(jTime, :),  ...
                        %obj.data{jChannel}.numMT(jTime, :), ...
                        %obj.data{jChannel}.mtRadialInt_raw(jTime,:),...
                        %obj.data{jChannel}.mtRadialInt(jTime,:),...
                        %obj.data{jChannel}.mtRadialVal(jTime,:)] = obj.analyzeMonopolar( mainFeature, jTime);

                case 'IMTBank'
                    %error('analyzeFeatureMicrotubule: InterphaseBank still in development')
                otherwise
                    error('analyzeFeatureMicrotubule: unknown mainFeature type')
            end


        end
        % }}}
        % analyzeSpindle {{{
        function spindleLength = analyzeSpindle( obj, mainFeature, jTime)
            % analyze a microtubule spindle            

            % Find spindle lengths
            spindleLength = norm( (mainFeature.featureList{1}.startPosition - mainFeature.featureList{1}.endPosition) .* obj.sizeVoxels );

        end
        % }}}
        % analyzeMonopolar {{{
        function [pole, numMT, radialInt_raw, radialInt, radialVal] = analyzeMonopolar( obj, mainFeature, jTime, data)
            % analyze a microtubule spindle            

            % Find monopolar pole position
            pole = mainFeature.featureList{1}.featureList{1}.position;

            % Find microtubule number
            numMT = mainFeature.featureList{1}.numFeatures - 1;

             %{
             {% Distance to Furthest Microtubule Tip
             {dst = 0;
             {for jmt = 1 : mainFeature.featureList{1}.numFeatures-1
             {    endDst = norm( mainFeature.featureList{1}.featureList{1+jmt}.startPosition - mainFeature.featureList{1}.featureList{1+jmt}.endPosition );
             {    if endDst > dst, dst = endDst; end 
             {end
             %}

            % get mt distribution away from pole
            if jTime == 50
                disp('here')
            end
            
            img = mat2gray( max(mainFeature.image, [], 3) );
            [radialInt_raw, radialVal] = Cell.phiIntegrate2D( img, pole(1:2), 30);
            % subtract background noise
            radialInt_raw = radialInt_raw - median( img( img > 0) );

            % normalize intensity distribution
            radialInt = radialInt_raw / max(radialInt_raw);

        end
        % }}}
        % graphChannelMicrotubule {{{
        function graphChannelMicrotubule( obj, jChannel)
            % graph analysis plots for the microtubule channel
            
            switch obj.cellType
                case 'Mitosis'
                    obj.graphSpindleLength( obj.data{jChannel}.spindleLength);

                case 'Interphase'

                case 'Monopolar'
                    obj.graphMicrotubuleNumber( obj.data{jChannel}.numMT);
                    obj.graphDistanceToTip( obj.data{jChannel}.tipDistance);

            end

        end
        % }}}
        % graphSpindleLength {{{
        function graphSpindleLength( obj, dat)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'spindle_length'); 
            plot( times, dat, 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Time (seconds)')
            ylabel( 'Length (microns)')
            set(gca, 'FontSize', 20)
            title('Spindle length vs time')
            hold off;
            drawnow 
            saveas( f, [obj.path, filesep, 'spindle_length.png'])
            close(f)

        end
        % }}}
        % graphMicrotubuleNumber {{{
        function graphMicrotubuleNumber( obj, dat)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'mt_number'); 
            plot( times, dat, 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Time (seconds)')
            ylabel( 'MT Number')
            set(gca, 'FontSize', 20)
            title('MT number vs time')
            hold off;
            drawnow 
            saveas( f, [obj.path, filesep, 'mt_number.png'])
            close(f)

        end
        % }}}
        % graphDistanceToTip{{{
        function graphDistanceToTip( obj, dat)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'mt_tip_distance'); 
            plot( times, dat, 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Time (seconds)')
            ylabel( 'Distance to Tip of Longest MT (voxels)')
            set(gca, 'FontSize', 20)
            title('Tip Distance vs time')
            hold off;
            drawnow 
            saveas( f, [obj.path, filesep, 'mt_tip_distance.png'])
            close(f)

        end
        % }}}
       % }}} 
       
        % Cut7 {{{
        % analyzeFeatureCut7 {{{
        function analyzeFeatureCut7( obj, mainFeature, jChannel, jTime)
            % do cut7 analysis
            
            switch obj.cellType
                case 'Mitosis'

                    [ obj.data{jChannel}.cut7SpindleInt(jTime,:), ...
                        obj.data{jChannel}.cut7SpindleInt_raw(jTime,:), ...
                        obj.data{jChannel}.cut7Range] = obj.analyzeCut7_Spindle( mainFeature, jTime);

                case 'Monopolar'

                    [ obj.data{jChannel}.cut7RadialInt(jTime,:), ...
                        obj.data{jChannel}.cut7RadialInt_raw(jTime,:), ...
                        obj.data{jChannel}.cut7RadialVal(jTime,:)] = obj.analyzeCut7_Monopolar( mainFeature, jTime);

                case 'Interphase'
                    error('analyzeFeatureMicrotubule: InterphaseBank still in development')
                otherwise
                    error('analyzeFeatureMicrotubule: unknown mainFeature type')
            end

        end
        % }}}
        % analyzeCut7_Spindle {{{
        function [spindleInt, spindleInt_raw, normRange] = analyzeCut7_Spindle( obj, mainFeature, jTime)
            
            cTime = obj.times( jTime);

            % get spindle start, spindle end from microtubule channel
            startPos = mainFeature.spindlePositionStart;
            endPos = mainFeature.spindlePositionEnd;

            % get cut7 frame 
            img = mainFeature.image;
            cut7amp = Cell.findAmplitudeAlongLine( max( img, [], 3), startPos(1:2), endPos(1:2));

            % normalize it from 0 to 1
            normRange = 0:0.02:1;

            % Store amplitude data normalized against spindle length
            spindleInt_raw = interp1( linspace(0,1,length( cut7amp) ), cut7amp, normRange);
            maxVal = max( spindleInt_raw); 
            spindleInt = spindleInt_raw / maxVal;

        end
        % }}}
        % analyzeCut7_Monopolar {{{
        function [radialInt, radialInt_raw, radialVal] = analyzeCut7_Monopolar( obj, mainFeature, jTime)
            
            % get pole of monopolar spindle
            pole = mainFeature.spindlePositionStart;

            % get cut7 frame  2D max projection
            imgZ = sum( mainFeature.image, 3);
            
            if jTime == 50
                disp('here')
            end
            
            % get mt sim image if it exists
%             if ~isempty( obj.simImageMT)
%                 imMT = obj.simImageMT( :,:,:,jTime);
%                 imMT = max( imMT, [], 3);
%                 imMT( imMT < 0.01*max(imMT(:)) ) = 0;
%                 imMT = logical( imMT);
%                 imgZ = imgZ .* imMT;
%             end

            % get cut7 distribution away from pole
            [radialInt_raw, radialVal] = Cell.phiIntegrate2D( imgZ, pole(1:2), 30);
            % subtract background noise
            radialInt_raw = radialInt_raw - median( imgZ( imgZ > 0) );

            % normalize intensity distribution
            radialInt = radialInt_raw / max(radialInt_raw);

        end
        % }}}
        % graphChannelCut7 {{{
        function graphChannelCut7( obj, jChannel)

            switch obj.cellType
                case 'Mitosis'

                    % graph cut7 intensity along the spindle averaged over time
                    obj.graphCut7MeanIntensity( obj.data{jChannel}.cut7SpindleInt, obj.data{jChannel}.cut7NormRange);
                    
                    % graph cut7 intensity along the spindle with 3rd dimension corresponding to time 
                    obj.graphCut7IntensityVsTime_Overlay( obj.data{jChannel}.cut7SpindleInt, obj.data{jChannel}.cut7NormRange);
                    obj.graphCut7IntensityVsTime_Surf( obj.data{jChannel}.cut7SpindleInt, obj.data{jChannel}.cut7NormRange);

                
                case 'Monopolar'
                    obj.graphCut7RadialIntensity_Monopolar( obj.data{jChannel}.cut7RadialInt, obj.data{jChannel}.cut7RadialVal);

            end

        end
        % }}}
        % graphCut7MeanIntensity{{{
        function graphCut7MeanIntensity( obj, ydat, xdat)

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_avg'); 
            plot( xdat, mean( ydat, 1), 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Distance along spindle ( normalized)')
            ylabel( 'Mean intensity')
            set(gca, 'FontSize', 20)
            title('Cut7 average distribution along spindle ')
            drawnow

            saveas( f, [obj.path, filesep, 'cut7_intensity_avg.png'])
            close(f)

        end
        % }}}
        % graphCut7IntensityVsTime_Overlay {{{
        function graphCut7IntensityVsTime_Overlay( obj, ydat, xdat)

            times = [ 1 : length(obj.times) ] .* mean(obj.timeStep);

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_overlay'); 
            hold on;
            for j = 1 : length(times)
                plot( xdat, ydat(j, :) );
            end
            hold off
            xlabel( 'Distance along spindle ( normalized)')
            ylabel( 'Mean intensity')
            set(gca, 'FontSize', 20)
            title('Cut7 distribution along spindle overlay')
            drawnow

            saveas( f, [obj.path, filesep, 'cut7_intensity_spindle_time_overlay.png'])
            close(f)

        end
        % }}}
        % graphCut7IntensityVsTime_Surf {{{
        function graphCut7IntensityVsTime_Surf( obj, ydat, xdat)

            times = [ 1 : length(obj.times) ] .* mean(obj.timeStep);
            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_surf'); 
            surf( xdat, times, ydat, 'LineStyle', ':', 'FaceColor', 'interp');
            xlabel( 'Distance along spindle ( normalized)')
            ylabel( 'Time')
            zlabel( 'Mean intensity')
            set(gca, 'FontSize', 20)
            title('Cut7 distribution along spindle with time')
            drawnow

            saveas( f, [obj.path, filesep, 'cut7_intensity_spindle_time_surf.png'])
            close(f)

        end
        % }}}
        % graphCut7RadialIntensity_Monopolar {{{
        function graphCut7RadialIntensity_Monopolar(obj, ydat, xdat)

            f = figure('NumberTitle', 'off', 'Name', 'cut7_radial_intensity_monopolar'); 
            plot( xdat, mean( ydat,1) , 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Distance from pole ( voxels)')
            ylabel( 'Total Intensity')
            set(gca, 'FontSize', 20)
            title('Cut7 intensity vs distance from pole')
            drawnow

            saveas( f, [obj.path, filesep, 'cut7_radial_intensity_monopolar.png'])
            close(f)

        end
        % }}}
        % }}}
        
        % findTimes {{{
        function times = findTimes( obj, channel)
            % scan the specific results path to get as list of times for which data exists 
                
            addpath( genpath( obj.path) )

            % get all folders names starting with C1 or C2 or C3 etc
            timeFolderFormat = ['C', num2str( channel) , '_*'];
            timeFolders = dir( [obj.path, filesep, timeFolderFormat]);
            timeFolders = {timeFolders.name};

            % Find start and end times from folder names
            times = erase( timeFolders, ['C' num2str( channel) '_T']);
            times = sort( str2double( times) );

        end
        % }}}

        % graphChannel {{{
        function graphChannel( obj, jChannel)
            % make analysis graphs for a single channel

            channel = obj.channels( jChannel);

            switch obj.features{jChannel}

                case 'Microtubule'

                    obj.graphChannelMicrotubule(jChannel )

                case 'Cut7'

                    obj.graphChannelCut7( jChannel)

            end

        end
        % }}}

        % makeMovie {{{
        function makeMovie( obj, jChannel)
            % Create frames for a movie
            
            global COUNTER
            COUNTER =1;
            scale_bars = 1;
            custom_xy = 1;
            for jTime = 1 : length( obj.times)
                
                % get feature
                feat = obj.data{jChannel}.features{jTime};
                feat.ID = COUNTER;
                feat.updateFeatureIDs();
                feat.updateFeatureMap();

                img = im2double(feat.image);
                im2 = max(img,[],3);
                
                % Get budding yeast mask
                if strcmp(obj.cellType, 'MitosisBud')
                    % Get a smoothed 2D image of the budding yeast cell. Smoothing allows
                    % values inside the cytoplasm to become similar.
                    im22 = imgaussfilt( max(img,[],3), 3);
                    % Threshold into 4 levels (bkg, cyto, mts, spindle)
                    tt = multithresh( im22, 3);
                    % Binarize the image, and fill all the holes
                    bw = imbinarize( im22, tt(1) );
                    se = strel('disk', 2,6);
                    bw = imerode( imdilate( bw, se), se);
                    bw = imfill(bw, 'holes');

                    % Smooth out the boundaries of this mask using convolution
                    windowSize = 11;
                    kernel = ones(windowSize) / windowSize ^ 2;
                    blurryImage = conv2(single(bw), kernel, 'same');
                    binaryImage = blurryImage > 0.5;

                    % Pick the biggest connected component
                    [labeledImage, numberOfBlobs] = bwlabel(binaryImage);
                    blobMeasurements = regionprops(labeledImage, 'area');
                    % Get all the areas
                    allAreas = [blobMeasurements.Area];
                    % Sort in order of largest to smallest.
                    [sortedAreas, sortIndexes] = sort(allAreas, 'descend');
                    % Extract the largest blob using ismember().
                    biggestBlob = ismember(labeledImage, sortIndexes(1));
                    % Convert from integer labeled image into binary (logical) image.
                    binaryImage = biggestBlob > 0;
                    imMask = binaryImage;
                    imMask = ones( size(im2));
                    im2 = im2.*imMask;
                end
                % Get xand y ranges for non-mask pixels
%                 xrange = find( any(im2,1)); yrange = find( any(im2,2));
                xrange = 1: size(im2,1); yrange = 1:size(im2,2);
                if custom_xy
                    if strcmp(obj.cellType, 'MitosisBud')
                        yrange = 1:size(im2,2); xrange = 1:123;
                    elseif strcmp(obj.cellType, 'Mitosis')
                        yrange = 38:113; xrange = 38:113;
                    elseif strcmp(obj.cellType, 'Monopolar')
                        yrange = 38:113; xrange = 38:113;
                    end
                end
                
                fs = 12;
                % make figure
                f = figure('visible', 'on'); 
                h = tight_subplot(1,3, 0.01, 0.3, 0.1);
                % Axes (1,1): Original
                set(f, 'currentaxes', h(1) );
                
                % Smooth, contrast stretch image for display
                J = im2; mask = im2; mask( mask >0) = 1;
                th = multithresh(im2(:), 3);
                J( J(:) > th(end) ) = th(end);
%                 J = J-median( J( mask(:)>0 ) );
%                 J( J(:) < 0) = 0;
%                 tols = [0.01, 0.9];
%                 pixes = sort( im2(im2(:) > 0) );
%                 tolval = pixes( round( tols*length(pixes)) );
%                 J = im2;
%                 J( J < tolval(1)) = min(im2(:));
%                 J( J > tolval(2)) = max(im2(:));
                J = imgaussfilt( J,1).*mask;
%                 J=im2;
%                 J = imgaussfilt( im2,1).*mask;
                imagesc( h(1), J ), 
                colormap( h(1), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); set( h(1), 'xtick', [], 'ytick', []);
                hold on;
                %title('Original Image');
                
                % Scale bar
                if scale_bars
                    PixelSize = obj.sizeVoxels(1); Scalebar_length = 1;
                    xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
                    % x_location and y_location are wherever you want your scale bar to appear.
                    line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
                end
                
                % Axes (2,1): Simulated Image
                set(f, 'currentaxes', h(2) );
                if strcmp(obj.cellType, 'MitosisBud')
                    imagesc( h(2), imMask.*max( feat.simulateAll( img, feat.ID) , [], 3) ),
                else
                    imagesc( h(2), max( feat.simulateAll( img, feat.ID) , [], 3) ), 
                end
                colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); set( h(2), 'xtick', [], 'ytick', []);
                %title('Fitted Image');
                
                % Axes (2,2): 3D features
                % Axes 1 for image background
                set(f, 'currentaxes', h(3) );hold on
                imagesc( h(3), J ), colormap(h(3), gray);hold on; h(3).Color = 'Black';
                axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); set( h(3), 'xtick', [], 'ytick', []);
                %title('3D Features');
                
                % Axes 2 for 3D features
                h_temp = axes('Position',get(h(3),'Position')); % make new axes
                set(h_temp, 'Color', 'none'); % i thought this may make the new axes background transparent, but it doesn't work 

                axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
                set( h_temp, 'xtick', [], 'ytick', []);
                set(f, 'currentaxes', h_temp );hold on
                colormap( h_temp, hsv);
                feat.displayFeature( h_temp,1);
                suptitle( sprintf('T=%d',obj.times(jTime)))
                
                % Axes 3 for colorbar
                h_temp2 = axes('Position',get(h(3),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
                axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
                set(f, 'currentaxes', h_temp2 );hold on
                colormap( h_temp2, hsv);

                tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
                tcks_label = ([1:2: size(img,3)]-1)*obj.sizeVoxels(3);
                cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
                    'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
                % cb.Label.String = 'Z(\mum)'; cb.Label.Rotation = 0;
                set(h_temp2, 'Position',get(h(3),'Position'))

                % Change cb height to be within 90% of axes
                pos0 = get(cb, 'Position'); pos1=pos0;
                cen = (pos0(2) + pos0(2)+pos0(4))/2;
                bottom = cen - 0.95*(cen-pos0(2)); top = cen + 0.95*(cen-pos0(2));
                pos1([2,4]) = [ bottom, top-bottom];
                pos1(3) = 0.015;
                set(cb, 'Position', pos1)
                
                % Set All font sizes
                for ii = 1:length(h)
                    set(h(ii),'FontSize',fs)
                end
                f.Color = 'white';
                obj.data{jChannel}.mov( jTime) = getframe(f);
                close(f)
                
                % TEMP CODE: 3D FIGURE
                time_enable = 109;
                if time_enable == obj.times(jTime)
                    stop_for_script=1;
                end

            end
            
            obj.writeMovie(jChannel);

        end
        % }}}

        % writeMovie {{{
        function writeMovie( obj, jChannel)
            % Write the frame to a movie file
            
            switch obj.features{jChannel}
                case 'Microtubule'
                    mname = 'mt_video';
                case 'Kinetochore'
                    mname = 'kc_video';
                case 'Cut7'
                    mname = 'cut7_video';
                otherwise
                    error('analyzeFrame: unknown feature')
            end

            writerObj = VideoWriter([obj.path, filesep, mname],'MPEG-4');
            writerObj.FrameRate = 10;
            writerObj.Quality = 100;

            % open the video writer
            open(writerObj);

            % write the frames to the video
            for frame = 1 : length( obj.data{jChannel}.mov)
                % convert the image to a frame
                writeVideo( writerObj, obj.data{jChannel}.mov(frame) );
            end

            % close the writer object
            close(writerObj);

        end
        % }}}
        
        function feat = linkForwardReverse( obj, feat1, feat2)
           
            % Accept the feature with the least residual
            sim1 = feat1.simulateAll( feat1.image, feat1.ID);
            sim2 = feat2.simulateAll( feat2.image, feat2.ID);
            r1 = sum( (feat1.image(:) - sim1(:) ).^2);
            r2 = sum( (feat2.image(:) - sim2(:) ).^2);
            if any(feat1.image(:) ~= feat2.image(:) )
                stophere = 1;
            end
            if abs(r1-r2)<0.02
                % pick the one with most tubulin
                switch feat1.type
                    case 'MonopolarAster'
                        L1 = 0; L2 = 0;
                        for j1 = 2 : feat1.featureList{1}.numFeatures
                            L1 = L1 + feat1.featureList{1}.featureList{j1}.length;
                        end
                        for j2 = 2 : feat2.featureList{1}.numFeatures
                            L2 = L2 + feat2.featureList{1}.featureList{j2}.length;
                        end
                        if L1 >= L2
                            feat = feat1.copyDeep();
                        else
                            feat = feat2.copyDeep();
                        end
                    case 'Spindle' || 'SpindleNew'
                        % get the one with the longest spindle
                        if feat1.featureList{1}.GetLength() >= feat2.featureList{1}.GetLength()
                            feat = feat1.copyDeep();
                        else
                            feat = feat2.copyDeep();
                        end
                    case 'IMTBank'
                        error('nada')
                end
            elseif r1 > r2
                feat = feat2.copyDeep();
            elseif r1 <= r2
                feat = feat1.copyDeep();
            end
            
%             feat = feat1.copyDeep();
%             switch feat1.type
%                 case 'MonopolarAster'
%                     n1 = feat1.featureList{1}.numFeatures; n2 = feat2.featureList{1}.numFeatures;
%                     
%                     f1 = @(j,x) feat1.featureList{1}.featureList{j}.(x);
%                     f2 = @(j,x) feat2.featureList{1}.featureList{j}.(x);
%                     % SPB averages
%                     feat.featureList{1}.featureList{1}.position = mean( ...
%                         [f1(1,'position'); f2(1,'position')],1);
%                     feat.featureList{1}.featureList{1}.amplitude = mean( ...
%                         [f1(1,'amplitude'), f2(1,'amplitude')]);
%                     feat.featureList{1}.featureList{1}.sigma = mean( ...
%                         [f1(1,'sigma'); f2(1,'sigma')],1);
%                     
%                     % Get params: length, theta, amplitude
%                     for jFeat = 1 : feat1.featureList{1}.numFeatures -1
%                         par1{jFeat}.amplitude = feat1.featureList{1}.featureList{1+jFeat}.amplitude;
%                         par1{jFeat}.length = feat1.featureList{1}.featureList{1+jFeat}.length;
%                         par1{jFeat}.theta = feat1.featureList{1}.featureList{1+jFeat}.theta;
%                         par1{jFeat}.theta(par1{jFeat}.theta < 0) = par1{jFeat}.theta(par1{jFeat}.theta < 0) + 2*pi;
%                     end
%                     for jFeat = 1 : feat2.featureList{1}.numFeatures -1
%                         par2{jFeat}.amplitude = feat2.featureList{1}.featureList{1+jFeat}.amplitude;
%                         par2{jFeat}.length = feat2.featureList{1}.featureList{1+jFeat}.length;
%                         par2{jFeat}.theta = feat2.featureList{1}.featureList{1+jFeat}.theta;
%                         par2{jFeat}.theta(par2{jFeat}.theta < 0) = par2{jFeat}.theta(par2{jFeat}.theta < 0) + 2*pi;
%                     end
%                     
%                     % find costs
%                     costs = zeros([n1-1 n2-1]);
%                     for j1 = 1:length(par1)
%                         for j2 = 1:length(par2)
%                             ad = par1{j1}.amplitude - par2{j2}.amplitude;
%                             ld = par1{j1}.length - par2{j2}.length;
%                             td = mod( abs(par1{j1}.theta - par2{j2}.theta),pi);
%                             costs(j1,j2) = (10*abs(ad))^2 + abs(ld)^2 + sum((5*abs(td)).^2);
%                         end
%                     end
%                     [pairs, s1, s2] = findLinks(costs);
%                     feat.featureList{1}.featureList = {feat.featureList{1}.featureList{1}};
%                     feat.featureList{1}.numFeatures = 1;
%                     for p = 1:size(pairs,1)
%                         ff1 = feat1.featureList{1}.featureList{1+pairs(p,1)};
%                         ff2 = feat2.featureList{1}.featureList{1+pairs(p,2)};
%                         ff = ff1.copy();
%                         ff.SetStartPosition( feat.featureList{1}.featureList{1}.position);
%                         ff.SetEndPosition( mean( [f1(1+pairs(p,1),'endPosition'); f2(1+pairs(p,2),'endPosition')],1) );
%                         ff.SetLength( mean( [f1(1+pairs(p,1),'length'); f2(1+pairs(p,2),'length')],1) );
%                         ff.SetTheta( mean( [f1(1+pairs(p,1),'theta'); f2(1+pairs(p,2),'theta')],1) );
%                         ff.amplitude = mean( [f1(1+pairs(p,1),'amplitude'); f2(1+pairs(p,2),'amplitude')],1);
%                         ff.sigma = mean( [f1(1+pairs(p,1),'sigma'); f2(1+pairs(p,2),'sigma')],1);
%                         feat.featureList{1}.addFeatureToList( ff);
%                     end
%                     for p = 1:length(s1)
%                         ff1 = feat1.featureList{1}.featureList{1+s1(p)};
%                         ff = ff1.copy();
%                         ff.SetStartPosition( feat.featureList{1}.featureList{1}.position);
%                         feat.featureList{1}.addFeatureToList( ff);
%                     end
%                     for p = 1:length(s2)
%                         ff1 = feat2.featureList{1}.featureList{1+s2(p)};
%                         ff = ff1.copy();
%                         ff.SetStartPosition( feat.featureList{1}.featureList{1}.position);
%                         feat.featureList{1}.addFeatureToList( ff);
%                     end
%                     
%                     % force inside z if any feature extends past z stack.
%                     if feat.dim == 3
%                         for jj = 2: feat.featureList{1}.numFeatures
%                             jf = feat.featureList{1}.featureList{jj};
% 
%                             outsideZ = ( jf.endPosition(3)<1 || jf.endPosition(3)>size(feat.image,3));
%                             while outsideZ
%                                 jf.SetLength( jf.length*0.99 );
%                                 outsideZ = ( jf.endPosition(3)<1 || jf.endPosition(3)>size(feat.image,3));
%                             end
%                         end         
%                     end
%                     
%                 case 'Spindle'
%                     error('not setup')
%                 case 'IMTBank'
%                     error('not setup')
%             end
            
            function [pairs, singles1, singles2] = findLinks(costs)
                % find links
                costAccept = 10;
                pairs = []; singles1 = []; singles2 = [];
                for j1 = 1 : size(costs,1)
                    [cmin, midx] = min( costs( j1,:) );
                    if cmin < costAccept && costs(midx, j1)
                        pairs = [ pairs; [j1,midx]];
                    else
                        singles1 = [ singles1(:), j1];
                    end
                end
                for j2 = 1 : size(costs,2)
                    %check if this has been matched
                    if ~any(pairs(:,2) == j2)
                        singles2 = [singles2 , j2];
                    end
                end
            end

        end
        
        % makeMovieFwdReverse {{{
        function makeMovieFwdReverse( obj, jChannel)
            % Create fwd reverse linking frames for movie
            
            if ~obj.fwdreverse
                return
            end
            
            global COUNTER
            COUNTER =1;
            
            for jTime = 1 : length( obj.times)
                
                % get forward feature
                featF = obj.data{jChannel}.featuresFwd{jTime};
                featF.ID = COUNTER;
                try
                featF.updateFeatureIDs();
                catch
                    stoph=1;
                end
                featF.updateFeatureMap();
                
                % get reverse feature
                featR = obj.data{jChannel}.featuresRev{jTime};
                featR.ID = COUNTER;
                featR.updateFeatureIDs();
                featR.updateFeatureMap();
                % set colors correctly
                AnalysisSingleCell.SetFeatureColor( featF, [0 1 0 0.5]); % Green
                AnalysisSingleCell.SetFeatureColor( featR, [0 0 1 0.5]); % Blue
                
                %get linked feature
                feat = obj.linkForwardReverse( featF, featR);
                feat.ID = COUNTER;
                feat.updateFeatureIDs();
                feat.updateFeatureMap();
                
                % Image
                img = im2double(featF.image);

                % make figure
                f = figure('visible', 'on'); 
                h = tight_subplot(2,2, [.05 .01]);
                
                % Axes (1,1): Original
                set(f, 'currentaxes', h(1) );
                imagesc( h(1), max(img , [], 3) ), 
                colormap( h(1), gray); axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(1), 'xtick', [], 'ytick', []);
                title('Original Image');
                
                % Axes (1,2): Original + 2D features
                set(f, 'currentaxes', h(2) );
                imagesc( h(2), max(img , [], 3) ), hold on
                colormap( h(2), gray); axis equal; axis ij; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(2), 'xtick', [], 'ytick', []);
                    
                featF.displayFeature( h(2));
                featR.displayFeature( h(2));
                title('2D Features');
                
                % Axes (2,1): Residual Comparison
                set(f, 'currentaxes', h(3) );
                simF = featF.simulateAll( img, featF.ID); simR = featR.simulateAll( img, featR.ID); simLinked = feat.simulateAll( img, feat.ID);
                resF = sum( (featF.image(:) - simF(:)).^2 ); resR = sum( (featR.image(:) - simR(:)).^2 ); resLinked = sum( (feat.image(:) - simLinked(:)).^2 );
                X = categorical({'Forward','Reverse','Linked'});
                Y = [ resF, resR, resLinked] - min([ resF, resR, resLinked]);
                bar(X,Y)
                title('Residual Comparison');
                
                % Axes (2,2): Linked 2D features
                set(f, 'currentaxes', h(4) );
                imagesc( h(4), max(img , [], 3) ), hold on
                colormap( h(4), gray); axis equal; axis ij; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(4), 'xtick', [], 'ytick', []);
                    
                feat.displayFeature( h(4));
                title('2D Linked Features');
                
                suptitle( sprintf('T=%d',obj.times(jTime)))
                obj.data{jChannel}.movFwdRev( jTime) = getframe(f);
                close(f)
            end
            
            obj.writeMovieFwdReverse(jChannel);

        end
        % }}}
        
        % writeMovieFwdReverse {{{
        function writeMovieFwdReverse( obj, jChannel)
            % Write the frame to a movie file
            
            switch obj.features{jChannel}
                case 'Microtubule'
                    mname = 'mt_video_fwd_reverse.avi';
                case 'Kinetochore'
                    mname = 'kc_video_fwd_reverse.avi';
                case 'Cut7'
                    mname = 'cut7_video_fwd_reverse.avi';
                otherwise
                    error('analyzeFrame: unknown feature')
            end

            writerObj = VideoWriter([obj.path, filesep, mname]);
            writerObj.FrameRate = 10;

            % open the video writer
            open(writerObj);

            % write the frames to the video
            for frame = 1 : length( obj.data{jChannel}.movFwdRev)
                % convert the image to a frame
                writeVideo( writerObj, obj.data{jChannel}.movFwdRev(frame) );
            end

            % close the writer object
            close(writerObj);

        end
        % }}}
        
        % makeMovie {{{
        function trackChannel( obj, jChannel)
            % Create frames for a movie
            
            % Get movie
            movieMat = zeros( size(obj.simImageMT,1), size(obj.simImageMT,2), size(obj.simImageMT,4) );
            for jTime = 1 : length( obj.times)
                feat = obj.data{jChannel}.features{jTime};
                movieMat(:,:,jTime) = max( im2double(feat.image),[],3);
            end
            
            % Get feature information
            switch obj.cellType
                case 'MitosisBud'
                    for jTime = 1 : length( obj.times)
                        feat = obj.data{jChannel}.features{jTime};
                        
                        % Create a structure with fields xCoord, yCoord,
                        % zCoord, amp, length, theta
                        % The coordinates will be of the endposition
                        % There is a single spindle with two asters inside.
                        xC = []; yC = []; zC = []; amp = []; el = []; phi = []; theta=[];
                        for ja = 2: feat.numFeatures
                            faster = feat.featureList{ja};
                            for jc = 2: faster.numFeatures
                                cc = faster.featureList{jc}.GetCoords();
                                xC = [xC; [cc(1,end), 0]];
                                yC = [yC; [cc(2,end), 0]];
                                zC = [zC; [cc(3,end), 0]];
                                amp = [amp; [faster.featureList{jc}.amplitude, 0]];
                                el = [el; [faster.featureList{jc}.L, 0]];
                                phi = [phi; [faster.featureList{jc}.thetaInit(1), 0]];
                                theta = [theta; [faster.featureList{jc}.thetaInit(2), 0]];
                            end
                        end
                        movieInfo(jTime).xCoords = xC;
                        movieInfo(jTime).yCoords = yC;
                        movieInfo(jTime).zCoords = zC;
                        movieInfo(jTime).amp = amp;
                        movieInfo(jTime).el = el;
                        movieInfo(jTime).theta = theta;
                        movieInfo(jTime).xCoord = el;
                        movieInfo(jTime).yCoord = phi;
                        movieInfo(jTime).zCoord = theta;
                            
                    end
                    tracksFinal = obj.trackMitosisBud(movieInfo);
                    % change tracksFinal to get coordinates instead of
                    % length, phi, theta
                    for idx = 1 : length(tracksFinal)
                        times = tracksFinal(idx).seqOfEvents(:,1);
                        for tt = times(1):times(2)
                            ii = tt+1-times(1);
                            track_number = tracksFinal(idx).tracksFeatIndxCG(ii);
                            if track_number == 0
                                continue
                            end
                            cm = [ movieInfo(tt).xCoords(track_number,:), ...
                                movieInfo(tt).yCoords(track_number,:),...
                                movieInfo(tt).zCoords(track_number,:),...
                                movieInfo(tt).amp(track_number,:),];
                            tracksFinal(idx).tracksCoordAmpCG( 8*ii-7: 8*ii) = cm;
                        end
                    end
                    plotTracks2D(tracksFinal, [], 2, [], 1, 1, zeros(140,123), [], 0, [], 1)
                    
                case 'Mitosis'
                    error('not set up')
                case 'Monopolar'
                    
                    startPt = zeros(3, length( obj.times) );
                    for jTime = 1 : length( obj.times)
                        feat = obj.data{jChannel}.features{jTime};
                        
                        % Create a structure with fields xCoord, yCoord,
                        % zCoord, amp, length, theta
                        % The coordinates will be of the endposition
                        % There is a single spindle with two asters inside.
                        xC = []; yC = []; zC = []; amp = []; el = []; phi = []; theta=[];
                        faster = feat.featureList{1};
                        startPt(:, jTime) = faster.featureList{1}.position;
                        for jc = 2: faster.numFeatures
                            cc = faster.featureList{jc}.endPosition;
                            xC = [xC; [cc(1), 0]];
                            yC = [yC; [cc(2), 0]];
                            zC = [zC; [cc(3), 0]];
                            amp = [amp; [faster.featureList{jc}.amplitude, 0]];
%                             el = [el; [faster.featureList{jc}.L, 0]];
%                             phi = [phi; [faster.featureList{jc}.thetaInit(1), 0]];
%                             theta = [theta; [faster.featureList{jc}.thetaInit(2), 0]];
                        end
                        movieInfo(jTime).xCoords = xC;
                        movieInfo(jTime).yCoords = yC;
                        movieInfo(jTime).zCoords = zC;
                        movieInfo(jTime).amp = amp;
%                         movieInfo(jTime).el = el;
%                         movieInfo(jTime).theta = theta;
                        movieInfo(jTime).xCoord = xC;
                        movieInfo(jTime).yCoord = yC;
                        movieInfo(jTime).zCoord = zC;
                            
                    end
                    tracksFinal = obj.trackMonopolar(movieInfo);
                    makeMovieTracking( tracksFinal, [],0,0,[],0,0,0,1,[],1,0,[],0,1,movieMat, [],1,0,'mp4_unix', startPt)
%                     plotTracks2D(tracksFinal, [], 2, [], 1, 1, zeros(150,150), [], 0, [], 1)
                    
            end
        end
        % }}}
        
        % temp_fig1 {{{
        function temp_fig1( obj, jChannel)
            % Create frames for a movie
            
            global COUNTER
            COUNTER =1;
            custom_xy = 1;
            tlist = 1 : length( obj.times);
            switch obj.cellType
                case 'Monopolar'
                    jTime = find(obj.times==58);
                case 'Mitosis'
                    jTime = find(obj.times==38);
                case 'MitosisBud'
                    jTime = find(obj.times==109);
            end
            
            f = figure('visible', 'on');
            h = tight_subplot(1,3, 0.005, 0.125,0.11);
            fs = 16;
            %set(f,'Position', get(groot,'Screensize'))
                
            % get feature
            feat = obj.data{jChannel}.features{ jTime };
            feat.ID = COUNTER;
            feat.updateFeatureIDs();
            feat.updateFeatureMap();
            img = im2double(feat.image);
            im2 = max(img,[],3);
                
            xrange = 1: size(feat.image,1); yrange = 1:size(feat.image,2);
            if custom_xy
                if strcmp(obj.cellType, 'MitosisBud')
                    yrange = 1:size(feat.image,2); xrange = 1:123;
                elseif strcmp(obj.cellType, 'Mitosis')
                    yrange = 30:105; xrange = 33:108;
                elseif strcmp(obj.cellType, 'Monopolar')
                    yrange = 38:113; xrange = 38:113;
                end
            end
                
            % Axes (1,1): Original
            set(f, 'currentaxes', h(1) );

            % Smooth, contrast stretch image for display
            J = im2;
            mask = im2; mask( mask(:) > 0) = 1;
            th = multithresh(J(J(:)>0),2);
            J( J(:) > th(end)) = th(end);
            J = imgaussfilt( J, 1).*mask;

            imagesc( h(1), J ), 
            colormap( h(1), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
            set( h(1), 'xtick', [], 'ytick', []);
            title(['Time = ' num2str( jTime)])
            hold on;
            %title('Original Image');
                
            % Scale bar
            PixelSize = obj.sizeVoxels(1); Scalebar_length = 1;
            xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
            % x_location and y_location are wherever you want your scale bar to appear.
            line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
                
            % Axes (2,1): Simulated Image
            set(f, 'currentaxes', h(2) );
            imagesc( h(2), max( feat.simulateAll( img, feat.ID) , [], 3) ), 
            colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
            set( h(2), 'xtick', [], 'ytick', []);
            %title('Fitted Image');

            % Axes (2,2): 3D features
            % Axes 1 for image background
            set(f, 'currentaxes', h(3) );hold on
            imagesc( h(3), J ), colormap(h(3), gray);hold on; h(3).Color = 'Black';
            axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
            set( h(3), 'xtick', [], 'ytick', []);
            %title('3D Features');
                
            % Axes 2 for 3D features
            h_temp = axes('Position',get(h(3),'Position')); % make new axes
            set(h_temp, 'Color', 'none');
            axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
            set( h_temp, 'xtick', [], 'ytick', []);
            set(f, 'currentaxes', h_temp );hold on
            colormap( h_temp, hsv);
            feat.displayFeature( h_temp,1);

            % Axes 3 for colorbar
            h_temp2 = axes('Position',get(h(3),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
            axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
            set(f, 'currentaxes', h_temp2 );hold on
            colormap( h_temp2, hsv);

            tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
            tcks_label = ([1:2: size(img,3)]-1)*obj.sizeVoxels(3);
            cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
                'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
            cb.Label.String = 'Z(\mum)'; cb.Label.Rotation = 90;
            set(h_temp2, 'Position',get(h(3),'Position'))

            % Change cb height to be within 90% of axes
            pos0 = get(cb, 'Position'); pos1=pos0;
            cen = (pos0(2) + pos0(2)+pos0(4))/2;
            bottom = cen - 0.9*(cen-pos0(2)); top = cen + 0.9*(cen-pos0(2));
            pos1([2,4]) = [ bottom, top-bottom];
            pos1(3) = 0.01;
            set(cb, 'Position', pos1)
            
            % Set All font sizes
            for ii = 1:length(h)
                set(h(ii),'FontSize',fs)
            end
            f.Color = 'white';
            
            % 3D FIGURE
            f2 = figure('visible', 'on');
            hold on; ax = gca; axis ij; set(ax, 'View',[-28,8]); f2.Color = 'white';
            set(ax,'zlim',[1 size(img,3)])
            feat.displayFeature3D( ax, size(img,3)); set(ax,'FontSize',fs)
            %xlim([60,90]); ylim( [60,100]); zlim([5,7]);
            tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
            tcks_label = 1:2: size(img,3);
            set(gca,'FontSize',20)
            pos = get(gca, 'Position'); 
            set(gca,'Position', [pos(1), pos(2), pos(3), pos(3)*size(img,1)/size(img,2)]);

            % Change z-limits if needed
            % set(gca,'zlim',[4.5,7.5]) % bipolar
            set(gca,'xlim',[20,80]) % byeast

            % XY planes
            zlim = get(gca,'zlim');
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');
            % [x, y] = meshgrid( xlim(1):1:xlim(2), ylim(1):1:ylim(2) ); % Generate x and y data
            [x, y] = meshgrid( xlim, ylim ); % Generate x and y data
            z = zeros( size(x) ); % Generate z data
            cols = 0.97*ones([size(z),3]);
            hold on
            for zz = zlim
                ztemp = z;
                ztemp(:) = zz;
                ss = surf(x,y,ztemp,cols, 'EdgeColor', [0.9,0.9,0.9]);
            end
            set(gca,'xlim', xlim, 'ylim', ylim);
            
            voxel_size = [0.1,0.1,0.5];
            ztick = get(gca,'ztick');
            ytick = get(gca,'ytick');
            xtick = get(gca,'xtick');
            ztl = (ztick - ztick(1))*voxel_size(3);
            ytl = (ytick - ytick(1))*voxel_size(2);
            xtl = (xtick - xtick(1))*voxel_size(1);

            set(gca,'xticklabels', xtl, 'yticklabels', ytl, 'zticklabels', ztl);
            xlabel = get(gca, 'xlabel'); set(xlabel, 'String', 'X (\mum)','FontWeight', 'normal')
            ylabel = get(gca, 'ylabel'); set(ylabel, 'String', 'Y (\mum)','FontWeight', 'normal')
            zlabel = get(gca, 'zlabel'); set(zlabel, 'String', 'Z (\mum)','FontWeight', 'normal')
            set(gcf,'Position', 1.1*get(gcf,'Position'))
        end
        
        % temp_fig3 {{{
        function temp_fig3( obj, jChannel)
            % Create frames for a movie
            
            global COUNTER
            COUNTER =1;
            custom_xy = 1;
            
            clearvars ttt
            ttt.monopolar = [11,12,13,14,15,16];
            ttt.bipolar = [60,61,62,63,64,65];
            ttt.byeast = [101,103,105,107,109,111];
            ttt.anaphase_elongation = [20,42,64,86,108,130];
            ttt = ttt.anaphase_elongation;
            ttt_real = [];
            for it = ttt
                ttt_real = [ttt_real, find(obj.times == it)];
            end
            ttt_old = ttt;
            ttt = ttt_real;
            
            f = figure('visible', 'on');
            h = tight_subplot(3,length(ttt), 0.001, 0.125,0.11);
            fs = 16;
            set(f,'Position', get(groot,'Screensize'))
            for jTime = 1 : length( ttt)
                
                % get feature
                feat = obj.data{jChannel}.features{ ttt( jTime) };
                feat.ID = COUNTER;
                feat.updateFeatureIDs();
                feat.updateFeatureMap();
                img = im2double(feat.image);
                im2 = max(img,[],3);
                
                xrange = 1: size(feat.image,1); yrange = 1:size(feat.image,2);
                if custom_xy
                    if strcmp(obj.cellType, 'MitosisBud')
                        yrange = 1:size(feat.image,2); xrange = 1:123;
                    elseif strcmp(obj.cellType, 'Mitosis')
                        yrange = 30:105; xrange = 33:108;
                        yrange = 30:105; xrange = 33:108; % Anaphase
                    elseif strcmp(obj.cellType, 'Monopolar')
                        yrange = 38:113; xrange = 38:113;
                    end
                end
                
                axn = jTime; 
                % Axes (1,1): Original
                set(f, 'currentaxes', h(axn) );
                
                % Smooth, contrast stretch image for display
                J = im2;
                mask = im2; mask( mask(:) > 0) = 1;
                th = multithresh(J(J(:)>0),2);
                J( J(:) > th(end)) = th(end);
                J = imgaussfilt( J, 1).*mask;

                imagesc( h(axn), J ), 
                colormap( h(axn), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
                set( h(axn), 'xtick', [], 'ytick', []);
                title(['Time = ' num2str( ttt_old(jTime))])
                hold on;
                %title('Original Image');
                
                % Scale bar
                if axn == 1
                    PixelSize = obj.sizeVoxels(1); Scalebar_length = 1;
                    xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
                    % x_location and y_location are wherever you want your scale bar to appear.
                    line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
                end
                
                % Axes (2,1): Simulated Image
                set(f, 'currentaxes', h(axn+1*length(ttt)) );
                imagesc( h(axn+1*length(ttt)), max( feat.simulateAll( img, feat.ID) , [], 3) ), 
                colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
                set( h(axn+1*length(ttt)), 'xtick', [], 'ytick', []);
                %title('Fitted Image');
                
                % Axes (2,2): 3D features
                % Axes 1 for image background
                set(f, 'currentaxes', h(axn+2*length(ttt)) );hold on
                imagesc( h(axn+2*length(ttt)), J ), colormap(h(axn+2*length(ttt)), gray);hold on; h(axn+2*length(ttt)).Color = 'Black';
                axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
                set( h(axn+2*length(ttt)), 'xtick', [], 'ytick', []);
                %title('3D Features');
                
                % Axes 2 for 3D features
                h_temp = axes('Position',get(h(axn+2*length(ttt)),'Position')); % make new axes
                set(h_temp, 'Color', 'none');
                axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
                set( h_temp, 'xtick', [], 'ytick', []);
                set(f, 'currentaxes', h_temp );hold on
                colormap( h_temp, hsv);
                feat.displayFeature( h_temp,1);

                if jTime == length(ttt)
                    % Axes 3 for colorbar
                    h_temp2 = axes('Position',get(h(axn+2*length(ttt)),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
                    axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
                    set(f, 'currentaxes', h_temp2 );hold on
                    colormap( h_temp2, hsv);

                    tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
                    tcks_label = ([1:2: size(img,3)]-1)*obj.sizeVoxels(3);
                    cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
                        'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
                    % cb.Label.String = 'Z(\mum)'; cb.Label.Rotation = 0;
                    set(h_temp2, 'Position',get(h(axn+2*length(ttt)),'Position'))

                    % Change cb height to be within 90% of axes
                    pos0 = get(cb, 'Position'); pos1=pos0;
                    cen = (pos0(2) + pos0(2)+pos0(4))/2;
                    bottom = cen - 1.01*(cen-pos0(2)); top = cen + 1.01*(cen-pos0(2));
                    pos1([2,4]) = [ bottom, top-bottom];
                    pos1(3) = 0.01;
                    set(cb, 'Position', pos1)
                end
            end
            
            % Set All font sizes
            for ii = 1:length(h)
                set(h(ii),'FontSize',fs)
            end
            f.Color = 'white';
                            
        end
        
        % trackMitosisBud {{{
        function tracksFinal = trackMitosisBud( obj, movieInfo)
           
            gapCloseParam.timeWindow = 5; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

            %optional input:
            gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            %function name
            costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

            %parameters

            parameters.linearMotion = 0; %use linear motion Kalman filter.
            parameters.minSearchRadius = 6; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            parameters.maxSearchRadius = 6; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

            parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

            parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            % parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

            %optional input
            parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

            costMatrices(1).parameters = parameters;
            clear parameters
            
            %function name
            costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

            %parameters

            %needed all the time
            parameters.linearMotion = 0; %use linear motion Kalman filter.

            parameters.minSearchRadius = 6; %minimum allowed search radius.
            parameters.maxSearchRadius = 6; %maximum allowed search radius.
            parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

            parameters.brownScaling = [0.25 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            % parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

            parameters.ampRatioLimit = []; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

            parameters.lenForClassify = 3; %minimum track segment length to classify it as linear or random.

            parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

            parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

            parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            % parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
            parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

            parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

            %optional; if not input, 1 will be used (i.e. no penalty)
            parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^(n-1)).

            %optional; to calculate MS search radius
            %if not input, MS search radius will be the same as gap closing search radius
            parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

            %NEW PARAMETER
            parameters.gapExcludeMS = 1; %flag to allow gaps to exclude merges and splits

            %NEW PARAMETER
            parameters.strategyBD = -1; %strategy to calculate birth and death cost

            costMatrices(2).parameters = parameters;
            clear parameters
            
            kalmanFunctions.reserveMem  = 'kalmanResMemLM';
            kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
            kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
            kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

            %% additional input

            %saveResults
            saveResults.dir = '~/Documents/Projects/ImageAnalysis/SingleCell'; %directory where to save input and output
            saveResults.filename = 'tracksTest.mat'; %name of file where input and output are saved
            % saveResults = 0; %don't save results

            %verbose state
            verbose = 1;

            %problem dimension
            probDim = 3;

            %% tracking function call

            [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
                costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
            
        end
        % }}}
        
        % trackMonopolar {{{
        function tracksFinal = trackMonopolar( obj, movieInfo)
           
            gapCloseParam.timeWindow = 2; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            gapCloseParam.minTrackLen = 3; %minimum length of track segments from linking to be used in gap closing.

            %optional input:
            gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            %function name
            costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

            %parameters

            parameters.linearMotion = 1; %use linear motion Kalman filter.
            parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            parameters.maxSearchRadius = 5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

            parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

            parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            % parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

            %optional input
            parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

            costMatrices(1).parameters = parameters;
            clear parameters
            
            %function name
            costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

            %parameters

            %needed all the time
            parameters.linearMotion = 1; %use linear motion Kalman filter.

            parameters.minSearchRadius = 6; %minimum allowed search radius.
            parameters.maxSearchRadius = 6; %maximum allowed search radius.
            parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

            parameters.brownScaling = [0.25 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            % parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

            parameters.ampRatioLimit = []; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

            parameters.lenForClassify = 3; %minimum track segment length to classify it as linear or random.

            parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

            parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

            parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            % parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
            parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

            parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

            %optional; if not input, 1 will be used (i.e. no penalty)
            parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^(n-1)).

            %optional; to calculate MS search radius
            %if not input, MS search radius will be the same as gap closing search radius
            parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

            %NEW PARAMETER
            parameters.gapExcludeMS = 1; %flag to allow gaps to exclude merges and splits

            %NEW PARAMETER
            parameters.strategyBD = -1; %strategy to calculate birth and death cost

            costMatrices(2).parameters = parameters;
            clear parameters
            
            kalmanFunctions.reserveMem  = 'kalmanResMemLM';
            kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
            kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
            kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

            %% additional input

            %saveResults
            saveResults.dir = '~/Documents/Projects/ImageAnalysis/SingleCell'; %directory where to save input and output
            saveResults.filename = 'tracksTest.mat'; %name of file where input and output are saved
            % saveResults = 0; %don't save results

            %verbose state
            verbose = 1;

            %problem dimension
            probDim = 3;

            %% tracking function call

            [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
                costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
            
        end
        % }}}
        
    end

    methods ( Static = true, Access = public )

        % AnalyzeSingle {{{
        function anaCell = AnalyzeSingle( resPathCell, params)
            % uses AnalysisBank to analyze results from a single movie
            
            % Make sure we have the correct paths
            %warning('off', 'MATLAB:rmpath:DirNotFound');
            %rmpath( genpath(pwd) );
            %warning('on', 'MATLAB:rmpath:DirNotFound');
            addpath( pwd);
            params = initAnalysisParams();
%             addpath( genpath( params.pathParent) );
            addpath( genpath( [pwd, filesep, 'classes']) );

            % if params not provided, find the params
            if nargin < 2
                %params = initAnalysisParams();
                fprintf('\n')
                fprintf('-----------------------------------------------------------------\n')
                fprintf('------------------- ANALYSIS SINGLE CELL ------------------------\n')
                fprintf('-----------------------------------------------------------------\n\n')
            end

            % check if resPathCell is a fullpath, if not then make it full
            [pth, ~, ~] = fileparts( resPathCell);
            if isempty( pth)
                resPathCell = fullfile( params.pathParent, resPathCell);
            end
            [~,cellName,~] = fileparts( resPathCell);
            fprintf( 'Cell = %s\n', cellName);
           
            params.Cell = load( [resPathCell, filesep, 'params.mat']);
            params.resPath = resPathCell;

            % Initialize analysis object
            anaCell = AnalysisSingleCell( resPathCell, params.Cell.params.cellInfo.type, params.channelsToAnalyze, params.channelFeatures, params.Cell.timeStep, params.Cell.sizeVoxels, params.flagMovie, params.flagGraph);

            % analyze the cell
            anaCell.Analyze();

        end
        % }}}
        
        % AnalyzeMulti {{{
        function anaCells = AnalyzeMulti( resRegExp)
            % uses AnalysisBank to analyze results from multiple movies
            % second argument can also be a cell array of multiple
            
            % Make sure we have the correct paths
            %warning('off', 'MATLAB:rmpath:DirNotFound');
            %rmpath( genpath(pwd) );
            %warning('on', 'MATLAB:rmpath:DirNotFound');
            %addpath( pwd);
            addpath( genpath( [pwd, filesep, 'classes']) );
            params = initAnalysisParams();
            addpath( params.pathParent);

            % Get all folder names in parent
            f = dir( params.pathParent);
            folds = {f.name};
            folds = folds( [f.isdir]);

            % Remove all folders starting with dot 
            folds = folds( cellfun( @isempty, regexp( folds, '^[.]') ) );

            % Ensure that we only have the results folders in the correct format (start with YYMMDD_HHMM)
            folds = folds( cellfun( @(x) ~isempty(x), regexp( folds, '^\d{6}_\d{4}') ) );

            % match resRegExp
            folds = folds( cellfun( @(x) ~isempty(x), regexp( folds, resRegExp) ) );
            numCells = length( folds);

            fprintf('\n')
            fprintf('-----------------------------------------------------------------\n')
            fprintf('-------------------- ANALYSIS MULTI CELL ------------------------\n')
            fprintf('-----------------------------------------------------------------\n\n')
            fprintf('    Total Cells : %d\n', numCells);

            % for each cell we'll run an AnalyzeSingle function
            for jCell = 1 : numCells
                addpath( genpath( folds{jCell} ) );
                anaCells{jCell} = AnalysisSingleCell.AnalyzeSingle( folds{ jCell}, params ); 
            end

            % Graph comparisons
            % AnalysisSingleCell.GraphMultiCell( anaCells)

        end
        % }}}
        
        % GraphingMultiCell {{{
        function GraphMultiCell( Cells)
            % Perform Multi Cell Graphing
            
            fprintf('Comparing analyzed cells...\n')
            nCells = length( Cells);

            for jf = 1: length( Cells{1}.features)
                switch Cells{1}.features{jf}
                    case 'Microtubule'
                        channel_mt = jf;
                    case 'Cut7'
                        channel_cut7 = jf;
                    case 'KC'
                        channel_kc = jf;
                end
            end

            % Compare Spindle Cells
            switch Cells{1}.cellType
                case 'Spindle'

                    % Compare spindle length
                    AnalysisSingleCell.GraphMulti_SpindleLength( Cells, channel_mt)

                    % Compare Cut 7 distribution along the spindle
                    AnalysisSingleCell.GraphMulti_Cut7AlongSpindle( Cells, channel_cut7)

                case 'Monopolar'

                    % Compare MT  distribution away from pole
                    AnalysisSingleCell.GraphMulti_MTDistributionAwayFromPoleMicrons( Cells, channel_mt)

                    % Compare Cut 7 distribution away from pole
                    AnalysisSingleCell.GraphMulti_Cut7DistributionAwayFromPoleMicrons( Cells, channel_cut7)

            end
            
        end
        % }}}

        % GraphMulti_SpindleLength {{{
        function GraphMulti_SpindleLength( Cells, channel)
            % Compare spindle lengths between cells

            fprintf('Comparing spindle lengths...\n')
            nCells = length( Cells);

            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'spindle_lengths_all'); 
            ax = axes;
            grid on; grid minor
            xlabel( 'Time (seconds)')
            ylabel( 'Length (microns)')
            set( ax, 'FontSize', 20)
            title('Spindle lengths vs time')
            hold on;

            % Process data for plotting
            % Spindle lengths
            for jCell = 1 : nCells
                spindleLengths{jCell} = Cells{jCell}.data{channel}.spindleLength;
            end
            % Times
            for jCell = 1 : nCells
                times{jCell} = Cells{jCell}.times - Cells{jCell}.times(1);
            end
            
            % plot mean and shaded error
            AnalysisSingleCell.plotAreaShadedError( spindleLengths, times, ax);

            for jCell = 1 : nCells
                saveas( f, [Cells{jCell}.path, filesep, 'spindle_lengths_all.png'])
            end
            [par, ~, ~] = fileparts( Cells{1}.path );
            saveas( f, [par, filesep, 'spindle_lengths_all.png'])
            close(f)

        end
        % }}}

        % GraphMulti_Cut7AlongSpindle {{{
        function GraphMulti_Cut7AlongSpindle( Cells, channel)
            % Compare Cut7 distibutions along the spindle

            fprintf('Comparing cut7 distirbution along the spindle...\n')
            nCells = length( Cells);

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_avg_all'); 
            ax = axes;
            grid on; grid minor
            xlabel( 'Distance along spindle ( normalized)')
            ylabel( 'Mean intensity')
            set( ax, 'FontSize', 20)
            title('Cut7 average distribution along spindle ')
            hold on;

            % Process data for plotting
            % cut7 intensity 
            for jCell = 1 : nCells
                int = mean( Cells{jCell}.data{channel}.cut7SpindleInt, 1);
                if mean( int(1: floor(end/3) ) ) < mean( int( ceil(2*end/3): end) )
                    int = flip( int);
                end
                cut7int{jCell} = int/max(int);
            end
            % Distance along spindle 
            for jCell = 1 : nCells
                dist{jCell} = Cells{jCell}.data{channel}.cut7Range;
            end
            
            % plot mean and shaded error
            AnalysisSingleCell.plotAreaShadedError( cut7int, dist, ax);

            for jCell = 1 : nCells
                saveas( f, [Cells{jCell}.path, filesep, 'cut7_intensity_spindle_avg_all.png'])
            end
            [par, ~, ~] = fileparts( Cells{1}.path );
            saveas( f, [par, filesep, 'cut7_intensity_spindle_avg_all.png'])
            close(f)



        end
        % }}}

        % GraphMulti_MTDistributionAwayFromPoleMicrons{{{
        function GraphMulti_MTDistributionAwayFromPoleMicrons( Cells, channel)
            % Compare Cut7 distibutions along the spindle

            fprintf('Comparing tubulin distribution away from pole...\n')
            nCells = length( Cells);

            f = figure('NumberTitle', 'off', 'Name', 'tub_intensity_monopolar_all_um'); 
            ax = axes;
            grid on; grid minor
            xlabel( 'Distance from pole [\mum]', 'interpreter', 'Tex')
            ylabel( 'Mean intensity')
            set( ax, 'FontSize', 20)
            title('Tub intensity vs distance from pole')
            hold on;

            % Process data for plotting
            % cut7 intensity 
            for jCell = 1 : nCells
                mtint{jCell} = Cells{jCell}.data{channel}.mtRadialInt_raw;
                mtint{jCell} = mean( mtint{jCell}, 1);
                %cut7int{jCell} = cut7int{jCell}/ max( cut7int{jCell});
            end

            % Distance from pole 
            for jCell = 1 : nCells
                dist{jCell} = Cells{jCell}.data{channel}.mtRadialVal;
                dist{jCell} = mean( dist{jCell}, 1);
                dist{jCell} = dist{jCell}.*Cells{jCell}.sizeVoxels(1);
            end
            
            % plot mean and shaded error
            AnalysisSingleCell.plotAreaShadedError( mtint, dist, ax);

            for jCell = 1 : nCells
                saveas( f, [Cells{jCell}.path, filesep, 'tub_intensity_monopolar_all_um.png'])
            end
            [par, ~, ~] = fileparts( Cells{1}.path );
            saveas( f, [par, filesep, 'tub_intensity_monopolar_all_um.png'])
            close(f)
            save([par,filesep,'tubdist.mat'], 'mtint', 'dist', '-v7.3')

        end
        % }}}
        
        % GraphMulti_Cut7DistributionAwayFromPoleMicrons{{{
        function GraphMulti_Cut7DistributionAwayFromPoleMicrons( Cells, channel)
            % Compare Cut7 distibutions along the spindle

            fprintf('Comparing cut7 distribution away from pole...\n')
            nCells = length( Cells);

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_monopolar_all_um'); 
            ax = axes;
            grid on; grid minor
            xlabel( 'Distance from pole [\mum]', 'interpreter', 'Tex')
            ylabel( 'Mean intensity')
            set( ax, 'FontSize', 20)
            title('Cut7 intensity vs distance from pole')
            hold on;

            % Process data for plotting
            % cut7 intensity 
            for jCell = 1 : nCells
                try
                    cut7int{jCell} = Cells{jCell}.data{channel}.cut7RadialInt_raw;
                    cut7int{jCell} = mean( cut7int{jCell}, 1);
                catch
                    disp('here')
                    cut7int{jCell} = Cells{jCell}.data{channel}.cut7RadialInt_raw;
                    cut7int{jCell} = mean( cut7int{jCell}, 1);
                end
                %cut7int{jCell} = cut7int{jCell}/ max( cut7int{jCell});
            end

            % Distance from pole J
            for jCell = 1 : nCells
                dist{jCell} = Cells{jCell}.data{channel}.cut7RadialVal;
                dist{jCell} = mean( dist{jCell}, 1);
                dist{jCell} = dist{jCell}.*Cells{jCell}.sizeVoxels(1);
            end
            
            % plot mean and shaded error
            AnalysisSingleCell.plotAreaShadedError( cut7int, dist, ax);

            for jCell = 1 : nCells
                saveas( f, [Cells{jCell}.path, filesep, 'cut7_intensity_monopolar_all_um.png'])
            end
            [par, ~, ~] = fileparts( Cells{1}.path );
            saveas( f, [par, filesep, 'cut7_intensity_monopolar_all_um.png'])
            close(f)
            save([par,filesep,'cut7dist.mat'], 'cut7int', 'dist', '-v7.3')

        end
        % }}}
        
        % plotAreaShadedError {{{
        function plotAreaShadedError( ydata, xdata, ax)
            % plot shaded area error given a set of unevenly sampled ordered data with varying number of observations
            % ydata : cell array of vectors corresponding to samples (each sample can be of a different length)
            % xdata (optional) : cell array of vectors corresponding to x values (each sample can be of a different length but must match its counterpart in ydata)
            % ax (optional) : axes to plot in

            colorError = [0.7 0.2 0.5];
            alphaError = 0.2;

            colorMean = [0.7 0.2 0.5];
            alphaMean = 1;
            meanWidth = 3;

            colorSample = [0.7 0.2 0];
            alphaSample = 0.6;
            sampleWidth = 1; 

            % Ensure presence of axes
            if nargin < 3
                f = figure('NumberTitle', 'off', 'Name', 'plotAreaShadedError');
                ax = axes;
            end

            % Ensure presence of x data
            nSamp = length( ydata);
            if nargin < 2
                xdata = cellfun( @(x) 1:length(x), ydata, 'UniformOutput', false);
            end

            % We will use linear interpolation to make our observations evenly and finely spaced out

            % find smallest x-step
            xstep = min( cellfun( @(x) min( diff(x) ), xdata) );
            
            % find max and min extent of observation
            xmax = max( cellfun( @(x) x(end), xdata) );
            xmin = min( cellfun( @(x) x(1), xdata) );

            % construct a new fine xdata vector
            nSteps = 2*ceil( (xmax-xmin)/xstep );
            x = linspace( xmin, xmax, nSteps);

            % find ydata interpolated values for the new x data vector 
            % Note: values outside the domain of the original xdata are set to NaN
            for jSamp = 1 : nSamp
                ydata_even(jSamp, :) = interp1( xdata{jSamp}, ydata{jSamp}, x, 'linear');
            end

            % find the mean of the evenly sampled ydata
            for jT = 1 : length( x) 
                xValues = ydata_even(:,jT);
                goodValues = xValues( ~isnan( xValues ) );
                ymean(jT) = mean( goodValues);
                yerr(jT) = std( goodValues);
            end
            
            % construct the extended vectors to define a region for shading
            x_vec = [x, fliplr(x)];
            y_vec = [ ymean+yerr, fliplr( ymean-yerr) ];

            % Shade the error region
            patch = fill(x_vec, y_vec, colorError);
            set(patch, 'Edgecolor', 'none');
            set(patch, 'FaceAlpha', alphaError );

            % Plot the samples
            for jSamp = 1 : nSamp 
                line( xdata{jSamp}, ydata{jSamp}, 'LineWidth', sampleWidth, 'Color', [colorSample alphaSample], 'LineStyle', ':') 
            end
            hold on;
            
            % Plot the mean value
            plot( x, ymean, 'Color', [ colorMean, alphaMean], ...
                'LineWidth', meanWidth);
            hold off;

        end
        % }}}
        
        function SetFeatureColor( feat, col)
           % Change color of basic features in feat to col
           
            switch feat.type
                case 'Spindle'
                    setNewColor( feat.featureList{1}, col);
                    for m = 1 : length(feat.featureList)-1
                        for j = 1 : feat.featureList{m+1}.numFeatures
                            setNewColor( feat.featureList{m+1}.featureList{j}, col);
                        end
                    end
                case 'MonopolarAster'
                    for j = 1 : feat.featureList{1}.numFeatures
                        setNewColor( feat.featureList{1}.featureList{j}, col);
                    end                    
                case 'IMTBank'
                    for j = 1 : length(feat.featureList)
                        setNewColor( feat.featureList{j}, col);
                    end 
                   
            end
           
            function feet = setNewColor(feet, cool)
                idx = find( strcmpi(feet.display, 'Color'));
                feet.display{idx+1} = cool;
            end
            
            
        end

    end

end
