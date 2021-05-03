classdef AnalysisSingleCell < handle
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
                fname = [obj.path, filesep,'dydata.mat'];
                
                disp( 'temp bypass of analysis file check')
                fprintf('   Analysis: channel %d with feature %s...\n', obj.channels( jChannel), obj.features{ jChannel}); 
                dat = obj.analyzeChannel( jChannel); save( fname, 'dat', '-v7.3');

                % graph channel analysis
                if obj.flagGraph
                    fprintf('   Graphing...\n')
                    obj.graphChannel( jChannel);
                end

                % Make movie
                if obj.flagMovie 
                    fprintf('   Directing movies...\n')
%                     obj.temp_fig1( jChannel);
%                     obj.temp_fig3( jChannel);
%                     obj.makeMovie( jChannel);
                    if obj.fwdreverse
%                         obj.makeMovieFwdReverse( jChannel);
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
            
            fix_amp_bug = 1;
            if fix_amp_bug
                fprintf('Fixing amplitude bug...\n')               
                obj.fix_amp_scaling_bug();      
            end
            
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
            eval( ['mainFeature = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );
            
            mainFeature.ID = COUNTER;
            mainFeature.updateFeatureIDs();
            mainFeature.updateFeatureMap();
            obj.data{jChannel}.features{jTime} = mainFeature;
            
            % Run Feature Specific analysis
%             obj.analyzeFeature( mainFeature, jChannel, jTime);
            
%             switch mainFeature.type
%                 case 'SpindleNew'
%                     for ja = 2 : 3
%                         for jc = 2 : length(mainFeature.featureList{ja}.featureList)
%                             mainFeature.featureList{ja}.featureList{jc}.display = {'Color', [1 0 0] , 'LineWidth', 2};
%                         end
%                     end
%                 case 'Spindle'
%                     for ja = 2 : 3
%                         for jc = 2 : length(mainFeature.featureList{ja}.featureList)
%                             mainFeature.featureList{ja}.featureList{jc}.display = {'Color', [0 1 0] , 'LineWidth', 2};
%                         end
%                     end
%                 case 'MonopolarAster'
%                     for jc = 2 : length(mainFeature.featureList{1}.featureList)
%                         mainFeature.featureList{1}.featureList{jc}.display = {'Color', [0 1 0] , 'LineWidth', 2};
%                     end
%             end

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
                   mask = BY_find_nuclear_mask( feat.image);
                   mask = max(mask,[],3);
                   im2 = im2.*mask;
                end
                %Get xand y ranges for non-mask pixels
                xrange = find( any(im2,1)); yrange = find( any(im2,2));
                %xrange = 1: size(im2,2); yrange = 1:size(im2,1);
                if custom_xy
%                     yrange = 45:105; xrange = 45:105;
%                     yrange = 38:113; xrange = 38:113;
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
                imagesc( h(1), im2 ), 
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
                    imagesc( h(2), max( feat.simulateAll( img, feat.ID) , [], 3) ),
                else
                    imagesc( h(2), max( feat.simulateAll( img, feat.ID) , [], 3) ), 
                end
                colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); set( h(2), 'xtick', [], 'ytick', []);
                %title('Fitted Image');
                
                % Axes (2,2): 3D features
                % Axes 1 for image background
                set(f, 'currentaxes', h(3) );hold on
                imagesc( h(3), im2 ), colormap(h(3), gray);hold on; h(3).Color = 'Black';
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
        
        % makeMovie {{{
        function trackChannel( obj, jChannel)
            % Create frames for a movie
            
            % Flag for tracking movie
            do_movie_tracking=0;
            
            % Get movie
            movieMat = zeros( size(obj.simImageMT,1), size(obj.simImageMT,2), size(obj.simImageMT,4) );
            for jTime = 1 : length( obj.times)
                feat = obj.data{jChannel}.features{jTime};
                movieMat(:,:,jTime) = max( im2double(feat.image),[],3);
            end
            
            % Get feature information
            switch obj.cellType
                case 'MitosisBud'
                    
                    % Track the SPBs. We will track curves at each SPB separately.
                    idx1 = zeros( 1,length( obj.times) ); idx1(1) = 1;
                    idx2 = idx1; idx2(1) = 2;
                    
                    % for each timestep, we will find a global minimal cost
                    % of SPB assignment. This cost will scale as dist^2.
                    for jTime = 2 : length( obj.times)
                        % get positions at previous and current time
                        pos01 = obj.data{jChannel}.features{jTime-1}.featureList{1}.startPosition;
                        pos02 = obj.data{jChannel}.features{jTime-1}.featureList{1}.endPosition;
                        pos11 = obj.data{jChannel}.features{jTime}.featureList{1}.startPosition;
                        pos12 = obj.data{jChannel}.features{jTime}.featureList{1}.endPosition;
                        
                        cost11 = ( norm( pos11-pos01) + norm( pos12-pos02) )^2;
                        cost12 = ( norm( pos12-pos01) + norm( pos11-pos02) )^2;
                        % fprintf('Cost11 = %.2f\nCost12 = %.2f\n',cost11,cost12)
                        if cost11 <= cost12
                            idx1(jTime) = idx1(jTime-1);
                            idx2(jTime) = idx2(jTime-1);
                        else
                            idx1(jTime) = 1+mod( idx1(jTime-1),2);
                            idx2(jTime) = 1+mod( idx2(jTime-1),2);
                        end
                    end
                    
                    % get asters at each SPB
                    feats1 = cell(1, length( obj.times));
                    feats2 = feats1;
                    spindles = feats1;
                    for jTime = 1 : length( obj.times)
                        feats1{jTime} = obj.data{jChannel}.features{jTime}.featureList{ 1+ idx1(jTime) };
                        feats2{jTime} = obj.data{jChannel}.features{jTime}.featureList{ 1+ idx2(jTime) };
                        spindles{jTime} = obj.data{jChannel}.features{jTime}.featureList{ 1};
                    end
                    
                    % Tracks curves at first SPB
                    trackBud1 = TrackCurves( movieMat, obj.times, obj.timeStep, obj.sizeVoxels);
                    trackBud1 = trackBud1.parseMainFeature( feats1 );
                    trackBud1 = trackBud1.trackUTRACK();
                    if ~isempty(trackBud1.tracksFinal)
                        trackBud1 = trackBud1.parseTracksFinal(feats1);
                    end
                    
                    % Tracks curves at second SPB
                    trackBud2 = TrackCurves( movieMat, obj.times, obj.timeStep, obj.sizeVoxels);
                    trackBud2 = trackBud2.parseMainFeature( feats2 );
                    trackBud2 = trackBud2.trackUTRACK();
                    if ~isempty(trackBud2.tracksFinal)
                        trackBud2 = trackBud2.parseTracksFinal( feats2);
                    end
                    
                    if ~isempty( trackBud1.features) && ~isempty( trackBud2.features)
                        dymgmt = DynamicFeatureMgmt( { trackBud1.features{:}, trackBud2.features{:} } );
                    elseif ~isempty( trackBud1.features) && isempty( trackBud2.features)
                        dymgmt = DynamicFeatureMgmt( trackBud1.features);
                    elseif isempty( trackBud1.features) && ~isempty( trackBud2.features)
                        dymgmt = DynamicFeatureMgmt( trackBud2.features);
                    else
                        disp('No trackable features found. Returning out of tracking call.')
                        return
                    end
                    dymgmt.saveMat( obj.path)
                    dymgmt.saveCSV( obj.path)
                    
                    feats = {spindles, trackBud1, trackBud2};
                    
                case 'Mitosis'
                    error('not set up')
                case 'Monopolar'
                    
                    % get asters at  SPB
                    spbs = cell(1, length( obj.times));
                    for jTime = 1 : length( obj.times)
                        spbs{jTime} = obj.data{jChannel}.features{jTime}.featureList{ 1}.featureList{ 1};
                    end
                    
                    trackMono = TrackLines(movieMat, obj.times, obj.timeStep, obj.sizeVoxels);
                    trackMono = trackMono.parseMainFeature( obj.data{jChannel}.features );
                    trackMono = trackMono.trackUTRACK();
                    if ~isempty(trackMono.tracksFinal)
                        trackMono = trackMono.parseTracksFinal( obj.data{jChannel}.features );
                        writematrix(trackMono.nfeats_per_frame,[obj.path, filesep,'nfeats_per_frame.csv'])
                        dymgmt = DynamicFeatureMgmt( trackMono.features);
                        dymgmt.saveMat( obj.path)
                    end
                    feats = {spbs, trackMono};
                    
            end
            
            if do_movie_tracking
                TrackFeatures.makeMovie( obj.cellType, movieMat, obj.times, feats, obj.path );
            end
                    
        end
        % }}}
        
        % temp_fig1 {{{
        function temp_fig1( obj, jChannel)
            % Create frames for a movie
            
            global COUNTER
            COUNTER =1;
            custom_xy = 1;
            invert = 1;
            
            tlist = 1 : length( obj.times);
            switch obj.cellType
                case 'Monopolar'
                    jTime = find(obj.times==58);
                case 'Mitosis'
                    jTime = find(obj.times==38);
                case 'MitosisBud'
                    jTime = find(obj.times==111);
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
                    yrange = 1:85; xrange = 1:85;
                elseif strcmp(obj.cellType, 'Mitosis')
                    yrange = 40:90; xrange = 45:95;
                elseif strcmp(obj.cellType, 'Monopolar')
                    yrange = 50:100; xrange = 50:100;
                end
            end
                
            % Axes (1,1): Original
            set(f, 'currentaxes', h(1) );

            % Smooth, contrast stretch image for display
            imgauss = imgaussfilt(im2,1);
            med = median(imgauss(im2(:)~=0));
            J = imadjust(imgauss, [med, 3*med],[]);

            if invert
                imagesc( h(1), imcomplement(J) ), 
            else
                imagesc( h(1), J )
            end
            pause(0.1)
            colormap( h(1), gray); axis equal; xlim( h(1), [xrange(1) xrange(end) ]); ylim( h(1), [yrange(1) yrange(end) ]); 
            set( h(1), 'xtick', [], 'ytick', []);
            title(['Time = ' num2str( jTime)])
            hold on;
            %title('Original Image');
                
            % Scale bar
            PixelSize = obj.sizeVoxels(1); Scalebar_length = 1;
            xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
            % x_location and y_location are wherever you want your scale bar to appear.
            if invert
                line([xstart, xend],[y0,y0], 'Color','k', 'LineWidth', 4)
            else
                line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
            end
            % Axes (2,1): Simulated Image
            set(f, 'currentaxes', h(2) );
            if invert
                imagesc( h(2), imcomplement(max( feat.simulateAll( img, feat.ID) , [], 3)) ), 
            else
                imagesc( h(2), max( feat.simulateAll( img, feat.ID) , [], 3) ),
            end
             
            colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
            set( h(2), 'xtick', [], 'ytick', []);
            %title('Fitted Image');

            % Axes (2,2): 3D features
            % Axes 1 for image background
            set(f, 'currentaxes', h(3) );hold on
            if invert
                imagesc( h(3), imcomplement(J) ), 
            else
                imagesc( h(3), J )
            end
            colormap(h(3), gray);hold on; h(3).Color = 'Black';
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
            if strcmp(obj.cellType, 'MitosisBud')
                view_angle = [-50,32];
                zlimits = [5,11];
            elseif strcmp(obj.cellType, 'Mitosis')
                view_angle = [-17,22];
                zlimits = [4,7];
            elseif strcmp(obj.cellType, 'Monopolar')
                view_angle = [-50,32];
                zlimits = [1,7];
            end
            f2 = figure('visible', 'on');
            hold on; ax = gca; axis ij;  f2.Color = 'white';
            set(ax, 'View',[10,10]);
            set(ax,'zlim',zlimits)
            feat.displayFeature3D( ax, size(img,3)); set(ax,'FontSize',fs)
            tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
            tcks_label = 1:2: size(img,3);
            set(gca,'FontSize',20)
%             pos = get(gca, 'Position'); 
%             set(gca,'Position', [pos(1), pos(2), pos(3), pos(3)*size(img,1)/size(img,2)]);
            set(gca,'Position', [0.2, 0.2, 0.6,0.6]);
            
            % Change z-limits if needed
            xlimits = [xrange(1), xrange(end)];
            ylimits = [yrange(1), yrange(end)];
            
            % XY planes
            zlims = get(gca,'zlim');
            [x, y] = meshgrid( xlimits, ylimits ); % Generate x and y data
            z = zeros( size(x) ); % Generate z data
            cols = 0.97*ones([size(z),3]);
            hold on
            for zz = zlims
                ztemp = z;
                ztemp(:) = zz;
                ss = surf(x,y,ztemp,cols, 'EdgeColor', [0.9,0.9,0.9]);
            end
            set(gca,'xlim', xlimits, 'ylim', ylimits);
            
            voxel_size = [0.1,0.1,0.5];
            ztick = get(gca,'ztick');
            ytick = get(gca,'ytick');
            xtick = get(gca,'xtick');
%             ztl = (ztick - ztick(1))*voxel_size(3);
%             ytl = (ytick - ytick(1))*voxel_size(2);
%             xtl = (xtick - xtick(1))*voxel_size(1);
            ztl = (ztick )*voxel_size(3);
            ytl = (ytick )*voxel_size(2);
            xtl = (xtick )*voxel_size(1);

            set(gca,'xticklabels', xtl, 'yticklabels', ytl, 'zticklabels', ztl);
            xlabel = get(gca, 'xlabel'); set(xlabel, 'String', 'X (\mum)','FontWeight', 'normal')
            ylabel = get(gca, 'ylabel'); set(ylabel, 'String', 'Y (\mum)','FontWeight', 'normal')
            zlabel = get(gca, 'zlabel'); set(zlabel, 'String', 'Z (\mum)','FontWeight', 'normal')
            set(ax, 'View',view_angle);
%             set(gcf,'Position', 1.0*get(gcf,'Position'))
        end
        
        % temp_fig3 {{{
        function temp_fig3( obj, jChannel)
            % Create frames for a movie
            
            global COUNTER
            COUNTER =1;
            custom_xy = 1;
            invert = 1;
            
            clearvars ttt
            ttt.monopolar = [11,12,13,14,15,16];
            ttt.bipolar = [60,61,62,63,64,65];
            ttt.byeast = [101,103,105,107,109,111];
            ttt.anaphase_elongation = [30,55,80,105,130,155];
            ttt = ttt.byeast;
            ttt_real = [];
            for it = ttt
                ttt_real = [ttt_real, find(obj.times == it)];
            end
            ttt_old = ttt;
            ttt = ttt_real;
            
            f = figure('visible', 'on');
            h = tight_subplot(3,length(ttt), [0.002,-0.025], 0.125,0.11);
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
                
                imgauss = imgaussfilt(im2,1);
                med = median(imgauss(im2(:)~=0));
                J = imadjust(imgauss, [med, 4*med],[]);
                %figure; imshowpair(im2,J,'montage');
                
                xrange = 1: size(feat.image,1); yrange = 1:size(feat.image,2);
                if custom_xy
                    if strcmp(obj.cellType, 'MitosisBud')
                        yrange = 1:85; xrange = 1:85;
                    elseif strcmp(obj.cellType, 'Mitosis')
                        yrange = 40:90; xrange = 45:95;
                        yrange = 50:105; xrange = 50:105; % Anaphase
                    elseif strcmp(obj.cellType, 'Monopolar')
                        yrange = 50:100; xrange = 50:100;
                    end
                end
                
                axn = jTime; 
                % Axes (1,1): Original
                set(f, 'currentaxes', h(axn) );
                if invert
                    imagesc( h(axn), imcomplement(J) ),
                else
                    imagesc( h(axn), J ),
                end
                colormap( h(axn), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
                set( h(axn), 'xtick', [], 'ytick', []);
                title({['Frame = ' num2str( ttt_old(jTime))],['Time = ',num2str( (ttt(jTime)-ttt(1))*obj.timeStep),' sec']})
                hold on;
                %title('Original Image');
                
                % Scale bar
                if axn == 1
                    PixelSize = obj.sizeVoxels(1); Scalebar_length = 1;
                    xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
                    % x_location and y_location are wherever you want your scale bar to appear.
                    if invert
                        line([xstart, xend],[y0,y0], 'Color','k', 'LineWidth', 4)
                    else
                        line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
                    end
                end
                
                % Axes (2,1): Simulated Image
                set(f, 'currentaxes', h(axn+1*length(ttt)) );
                if invert
                    imagesc( h(axn+1*length(ttt)), imcomplement( max( feat.simulateAll( img, feat.ID) , [], 3)) ), 
                else
                    imagesc( h(axn+1*length(ttt)), max( feat.simulateAll( img, feat.ID) , [], 3) ), 
                end
                colormap gray; axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
                set( h(axn+1*length(ttt)), 'xtick', [], 'ytick', []);
                %title('Fitted Image');
                
                % Axes (2,2): 3D features
                % Axes 1 for image background
                set(f, 'currentaxes', h(axn+2*length(ttt)) );hold on
                if invert
                    imagesc( h(axn+2*length(ttt)), imcomplement(J) )
                else
                    imagesc( h(axn+2*length(ttt)), J )
                end
                colormap(h(axn+2*length(ttt)), gray);hold on; h(axn+2*length(ttt)).Color = 'Black';
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
        
        % fix_amp_scaling_bug 
        function obj = fix_amp_scaling_bug(obj)
            
            function new = scale_this_amp(old, oldMax,oldMin, newMax,newMin)
                new = (old-oldMin)*(newMax-newMin)/(oldMax-oldMin) + newMin;
            end
            function new = scale_this_amp_err(old, oldMax,oldMin, newMax,newMin)
                new = old*(newMax-newMin)/(oldMax-oldMin);
            end
            % Define path to search for amplitude file
            path_search = '~/Documents/Projects/ImageAnalysis/FY Datasets/Paper/Monopolar';
            
            % current cell name, try to find it.
            c_name = obj.folderName(14:end);
            afil = dir([path_search,filesep,'*/',c_name,'_amp.mat']);

            % load amplitude file if found
            if isempty(afil)
                fprintf('amp file not found\n')
            else
                adata = load([afil.folder,filesep,afil.name]);
                % Determine amplitude scaling
                % We have the true max/min as has been loaded
                % We will measure the current max/min from the images in
                % features
                nT = length(obj.times);
                maxAmp = adata.maxAmp;
                minAmp = adata.minAmp;
                maxAmp_c = zeros( size(maxAmp));
                minAmp_c = maxAmp_c;
                for jt = 1 : nT
                    img = obj.data{1}.features{jt}.image;
                    msk = obj.data{1}.features{jt}.mask;
                    vals = img( find(msk(:)) );
                    maxAmp_c(jt) = max( vals);
                    minAmp_c(jt) = min( vals);
                end
                % Data is for uint16, so we can say that:
                % 1 <-> 2^16 -1 = 65535
                % 0 <-> 0
                new_upper = maxAmp/65535;
                new_lower = minAmp/65535;
                for jt = 1 : nT
                    
                    imOld = obj.data{1}.features{jt}.image;
                    % scale image
                    imNew = obj.data{1}.features{jt}.mask.* scale_this_amp(imOld, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                    obj.data{1}.features{jt}.image = imNew;
                end
                % Scale all amplitudes, and backgrounds                
                for jt = 1 : nT
                    
                    % Main feature at current time
                    mf = obj.data{1}.features{jt};
                    
                    if isprop( mf, 'background')
                        mf.background = scale_this_amp(mf.background, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                    end
                        
                    for j1 = 1 : mf.numFeatures
                        
                        f1 = mf.featureList{j1};
                        
                        % Change amplitude, background
                        if isprop( f1, 'amplitude')
                            f1.amplitude = scale_this_amp(f1.amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                            f1.err_amplitude = scale_this_amp_err(f1.err_amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                        end
                        if isprop( f1, 'background')
                            f1.background = scale_this_amp(f1.background, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                        end
                        
                        % Should we go deeper into object
                        if ~isprop( f1, 'numFeatures')
                            continue;
                        end
                        
                        % Go deeper
                        for j2 = 1 : f1.numFeatures
                        
                            f2 = f1.featureList{j2};
                            
                            % Change amplitude, background
                            if isprop( f2, 'amplitude')
                                f2.amplitude = scale_this_amp(f2.amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                                f2.err_amplitude = scale_this_amp_err(f2.err_amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                            end
                            if isprop( f2, 'background')
                                f2.background = scale_this_amp(f2.background, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                            end

                            % Should we go deeper into object
                            if ~isprop( f2, 'numFeatures')
                                continue;
                            end
                            
                            % Go deeper
                            for j3 = 1 : f2.numFeatures

                                f3 = f2.featureList{j3};

                                % Change amplitude, background
                                if isprop( f3, 'amplitude')
                                    f3.amplitude = scale_this_amp(f3.amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                                    f3.err_amplitude = scale_this_amp_err(f3.err_amplitude, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                                end
                                if isprop( f3, 'background')
                                    f3.background = scale_this_amp(f3.background, maxAmp_c(jt),minAmp_c(jt), new_upper(jt), new_lower(jt) );
                                end
                            end
                        end 
                        
                    end
                end
                
                
            end
        end
        
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
                if exist([params.pathParent, filesep, folds{jCell}, filesep, 'nfeats_per_frame.csv'], 'file') ~=2
                    anaCells{jCell} = AnalysisSingleCell.AnalyzeSingle( folds{ jCell}, params ); 
                end
            end

            % Graph comparisons
            % AnalysisSingleCell.GraphMultiCell( anaCells)

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
