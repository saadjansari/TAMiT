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
    end

    methods (Access = public )

        % AnalysisSingleCell {{{
        function obj = AnalysisSingleCell( pathResCell, cellType, channels, features, timeStep, sizeVoxels, flagMovie, flagGraph)
            % Contruct the Analysis Object
            
            obj.path = pathResCell;
            if ~strcmp(cellType, 'Mitosis') && ~strcmp( cellType, 'Interphase') && ~strcmp( cellType, 'Monopolar')
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
                    obj.makeMovie( jChannel);
                end

            end

        end
        % }}}

        % analyzeChannel{{{
        function data = analyzeChannel( obj, jChannel)
            % Analyze this frame
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

            % Current channel 
            channel = obj.channels( jChannel);
            time = obj.times( jTime);

            % Load the data file if available
            timeDataFile = ['C' num2str( channel) '_T' num2str( time) '_final.mat']; 
            timeData = load( timeDataFile);
            obj.goodFrame = jTime;

            % Initialize main feature from loaded data
            eval( ['mainFeature = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );

            % Run Feature Specific analysis
            obj.analyzeFeature( mainFeature, jChannel, jTime);

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
            
            for jTime = 1 : length( obj.times)
                
                % get feature
                feat = obj.data{jChannel}.features{jTime};
                img = im2double(feat.image);

                % make figure
                f = figure('visible', 'on'); 
                h = tight_subplot(1,2);
                set(f, 'currentaxes', h(1) );
                imagesc( h(1), max(img , [], 3) ), 
                colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(1), 'xtick', [], 'ytick', []);
                title(sprintf('ref: T = %d', obj.times(jTime) ) );
                set(f, 'currentaxes', h(2) );
                imagesc( h(2), max(img , [], 3) ), hold on
                colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(2), 'xtick', [], 'ytick', []);
                feat.displayFeature( h(2) );
                title(sprintf('feature: T = %d', obj.times(jTime) ) );
                obj.data{jChannel}.mov( jTime) = getframe(f);
                close(f)

            end
            
            obj.writeMovie(jChannel);

        end
        % }}}

        % writeMovie {{{
        function writeMovie( obj, jChannel)
            % Write the frame to a movie file
            
            switch obj.features{jChannel}
                case 'Microtubule'
                    mname = 'mt_video.avi';
                case 'Kinetochore'
                    mname = 'kc_video.avi';
                case 'Cut7'
                    mname = 'cut7_video.avi';
                otherwise
                    error('analyzeFrame: unknown feature')
            end

            writerObj = VideoWriter([obj.path, filesep, mname]);
            writerObj.FrameRate = 10;

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
            AnalysisSingleCell.GraphMultiCell( anaCells)

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

    end

end
