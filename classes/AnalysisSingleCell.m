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
        data
        folderName
        goodFrame = 1
    end

    methods (Access = public )

        % AnalysisSingleCell {{{
        function obj = AnalysisSingleCell( pathResCell, cellType, channels, features, timeStep, sizeVoxels)
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

        end
        % }}}
        
        % Analyze {{{
        function Analyze( obj)
            % Analyze the object

            fprintf('Analyzing cell results for cell : %s\n', obj.folderName);

            % Analyze each channel
            for jChannel = 1: length( obj.channels)

                channel = obj.channels( jChannel);
                feature = obj.features{ jChannel};
                fprintf('   Analyzing channel %d with feature %s...\n', channel, feature); 

                % analyze a single channel
                obj.analyzeChannel( jChannel);

                % graph channel analysis
                obj.graphChannel( jChannel);

                % Make movie
                obj.makeMovie( jChannel);

            end

        end
        % }}}

        % analyzeChannel{{{
        function analyzeChannel( obj, jChannel)

            % Current channel 
            channel = obj.channels( jChannel);

            % find times for which analysis will be performed
            obj.times = obj.findTimes( channel);

            for jFrame = 1 : length( obj.times)

                % Current time frame
                jTime = obj.times( jFrame);
                fprintf('      Analyzing time %d...\n', jTime); 

                % Analyze this frame
                obj.analyzeFrame( jChannel, jFrame);

            end
            
        end
        % }}}

        % analyzeFrame {{{
        function analyzeFrame( obj, jChannel, jTime)
            % Analyze the time given in cTime

            % Current channel 
            channel = obj.channels( jChannel);
            time = obj.times( jTime);

            % Load the data file if available
            try
                timeDataFile = ['C' num2str( channel) '_T' num2str( time) '_final.mat']; 
                timeData = load( timeDataFile);
                obj.goodFrame = jTime;
            catch
                % If data file does not exist, use the last good frame available
                disp('oops')
                time = obj.times( obj.goodFrame);
                timeDataFile = ['C' num2str( channel) '_T' num2str( time) '_final.mat']; 
                timeData = load( timeDataFile);
            end

            % Initialize main feature from loaded data
            eval( ['mainFeature = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );

            % Run Feature Specific analysis
            obj.analyzeFeature( mainFeature, jChannel, jTime);

        end
        % }}}

        % analyzeFeature {{{
        function analyzeFeature( obj, mainFeature, jChannel, jTime)
            % Run analysis on the main feature by deciding what type of feature exists

            switch obj.features{jChannel}
                case 'Microtubule'
                    obj.analyzeFeatureMicrotubule( mainFeature, jTime);
                case 'Kinetochore'
                    obj.analyzeFeatureKinetochore( mainFeature, jTime);
                case 'Cut7'
                    obj.analyzeFeatureCut7( mainFeature, jTime);
                otherwise
                    error('analyzeFrame: unknown feature')
            end

        end
        % }}}
        
        % Microtubules {{{
        % analyzeFeatureMicrotubule {{{
        function analyzeFeatureMicrotubule( obj, mainFeature, jTime)
            % Analyze some type of microtubule feature (e.g. Spindle, Monopolar, InterphaseBank)
            
            switch mainFeature.type
                case 'Spindle'
                    % get spindle lengths
                    obj.data.spindle.mt = obj.analyzeSpindle( mainFeature, jTime);
                case 'MonopolarAster'
                    obj.data.monopolar.mt = obj.analyzeMonopolar( mainFeature, jTime);
                case 'InterphaseBank'
                    error('analyzeFeatureMicrotubule: InterphaseBank still in development')
                otherwise
                    error('analyzeFeatureMicrotubule: unknown mainFeature type')
            end

            % Capture frame for movie
            img = mainFeature.image;
            f = figure('visible', 'off'); 
            h = tight_subplot(1,2);
            set(f, 'currentaxes', h(1) );
            imagesc( h(1), max(img , [], 3) ), 
            colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(1), 'xtick', [], 'ytick', []);
            title(sprintf('ref: T = %d', obj.times(jTime) ) );
            set(f, 'currentaxes', h(2) );
            imagesc( h(2), max(img , [], 3) ), 
            colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(2), 'xtick', [], 'ytick', []);
            mainFeature.displayFeature( h(2) );
            title(sprintf('feature: T = %d', obj.times(jTime) ) );
            obj.data.movSpindle( jTime) = getframe( f);
            close(f)

        end
        % }}}
        % analyzeSpindle {{{
        function data = analyzeSpindle( obj, mainFeature, jTime)
            % analyze a microtubule spindle            

            % Find spindle lengths
            try; data = obj.data.spindle.mt; end
            data.spindleLength(jTime) = norm( (mainFeature.featureList{1}.startPosition - mainFeature.featureList{1}.endPosition) .* obj.sizeVoxels );

            % Any other things?

        end
        % }}}
        % analyzeMonopolar {{{
        function data = analyzeMonopolar( obj, mainFeature, jTime, data)
            % analyze a microtubule spindle            

            % Find monopolar pole position
            try; data = obj.data.monopolar.mt; end
            data.pole(jTime, :) = mainFeature.featureList{1}.featureList{1}.position;

            % Find microtubule number
            data.numMT(jTime) = mainFeature.featureList{1}.numFeatures - 1;

            % Distance to Furthest Microtubule Tip
            dst = 0;
            for jmt = 1 : mainFeature.featureList{1}.numFeatures-1
                endDst = norm( mainFeature.featureList{1}.featureList{1+jmt}.startPosition - mainFeature.featureList{1}.featureList{1+jmt}.endPosition );
                if endDst > dst, dst = endDst; end 
            end
            data.distEndTip(jTime) = dst;

        end
        % }}}
        % graphChannelMicrotubule {{{
        function graphChannelMicrotubule( obj)
            % graph analysis plots for the microtubule channel
            
            switch obj.cellType
                case 'Mitosis'
                    obj.graphSpindleLength();

                case 'Interphase'

                case 'Monopolar'
                    obj.graphMicrotubuleNumber();
                    obj.graphDistanceToTip();

            end

        end
        % }}}
        % graphSpindleLength {{{
        function graphSpindleLength( obj)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'spindle_length'); 
            plot( times, obj.data.spindle.mt.spindleLength, 'LineWidth', 2, 'Color', 'm') 
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
        function graphMicrotubuleNumber( obj)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'mt_number'); 
            plot( times, obj.data.monopolar.mt.numMT, 'LineWidth', 2, 'Color', 'm') 
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
        function graphDistanceToTip( obj)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'distance_mt_tip'); 
            plot( times, obj.data.monopolar.mt.distEndTip, 'LineWidth', 2, 'Color', 'm') 
            grid on
            grid minor
            xlabel( 'Time (seconds)')
            ylabel( 'Distance to Tip of Longest MT (voxels)')
            set(gca, 'FontSize', 20)
            title('Tip Distance vs time')
            hold off;
            drawnow 
            saveas( f, [obj.path, filesep, 'distance_mt_tip.png'])
            close(f)

        end
        % }}}
       % }}} 
       
        % Cut7 {{{
        % analyzeFeatureCut7 {{{
        function analyzeFeatureCut7( obj, mainFeature, jTime)
            % do cut7 analysis
            
            switch obj.cellType
                case 'Mitosis'
                    obj.data.spindle.cut7 = obj.analyzeCut7_Spindle( mainFeature, jTime);
                case 'Monopolar'
                    obj.data.monopolar.cut7 = obj.analyzeCut7_Monopolar( mainFeature, jTime);
                case 'Interphase'
                    error('analyzeFeatureMicrotubule: InterphaseBank still in development')
                otherwise
                    error('analyzeFeatureMicrotubule: unknown mainFeature type')
            end

            % Capture frame for movie
            img = mainFeature.image;
            startPos = mainFeature.spindlePositionStart;
            endPos = mainFeature.spindlePositionEnd;
            f = figure('visible', 'off'); 
            h = tight_subplot(1,2);
            set(f, 'currentaxes', h(1) );
            imagesc( h(1), max(img , [], 3) ), 
            colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(1), 'xtick', [], 'ytick', []);
            title(sprintf('ref: T = %d', obj.times(jTime) ) );
            set(f, 'currentaxes', h(2) );
            imagesc( h(2), max(img , [], 3) ), 
            colormap gray; axis equal; xlim( [1 size(img, 2) ]); ylim( [1 size(img, 1) ]); set( h(2), 'xtick', [], 'ytick', []);
            line( [startPos(1), endPos(1)], [startPos(2), endPos(2)], 'Color', [1 0.6 0.6 0.5], 'LineStyle', '-', 'LineWidth', 2)
            title(sprintf('feature: T = %d', obj.times(jTime) ) );
            obj.data.movCut7( jTime) = getframe( f);
            close(f)

        end
        % }}}
        % analyzeCut7_Spindle {{{
        function data = analyzeCut7_Spindle( obj, mainFeature, jTime)
            
            cTime = obj.times( jTime);

            % get spindle start, spindle end from microtubule channel
            startPos = mainFeature.spindlePositionStart;
            endPos = mainFeature.spindlePositionEnd;

            % get cut7 frame 
            img = mainFeature.image;
            cut7amp = Cell.findAmplitudeAlongLine( max( img, [], 3), startPos(1:2), endPos(1:2));

            % normalize it from 0 to 1
            try; data = obj.data.spindle.cut7; end
            normRange = 0:0.02:1;
            if ~isfield( obj.data, 'cut7Range');
                data.range = normRange;
            end

            % Store amplitude data normalized against spindle length
            data.spindleInt(jTime, :) = interp1( linspace(0,1,length( cut7amp) ), cut7amp, normRange);

        end
        % }}}
        % analyzeCut7_Monopolar {{{
        function data = analyzeCut7_Monopolar( obj, mainFeature, jTime)
            
            % get pole of monopolar spindle
            pole = mainFeature.spindlePositionStart;

            % get cut7 frame 
            img = mainFeature.image;
            % sum image in z direction
            imgSum = sum( img, 3);

            try; data = obj.data.monopolar.cut7; end
            % get cut7 distribution away from pole
            [data.radialInt_raw(jTime, :), data.radialVal(jTime, :)] = Cell.phiIntegrate2D( imgSum, pole(1:2) );
            data.radialInt(jTime, :) = data.radialInt_raw(jTime, :) / max( data.radialInt_raw(jTime, :) );

        end
        % }}}
        % graphChannelCut7 {{{
        function graphChannelCut7( obj)

            switch obj.cellType
                case 'Mitosis'
                    % normalize cut7 amplitude data over time
                    maxNorm = max( vecnorm( obj.data.spindle.cut7.spindleInt, 2, 2) );
                    maxVal = max( max( obj.data.spindle.cut7.spindleInt) );
                    obj.data.spindle.cut7.spindleInt_raw = obj.data.spindle.cut7.spindleInt;
                    obj.data.spindle.cut7.spindleInt = obj.data.spindle.cut7.spindleInt / maxVal;

                    % graph cut7 intensity along the spindle averaged over time
                    obj.graphCut7MeanIntensity;
                    
                    % graph cut7 intensity along the spindle with 3rd dimension corresponding to time 
                    obj.graphCut7IntensityVsTime_Overlay;
                    obj.graphCut7IntensityVsTime_Surf;
                
                case 'Monopolar'
                    
                    obj.graphCut7RadialIntensity_Monopolar;

            end

        end
        % }}}
        % graphCut7IntensityVsTime_Overlay {{{
        function graphCut7IntensityVsTime_Overlay( obj)

            times = [ 1 : length(obj.times) ] .* mean(obj.timeStep);
            distSpindle = obj.data.spindle.cut7.range;

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_overlay'); 
            hold on;
            for j = 1 : length(times)
                plot( distSpindle, obj.data.spindle.cut7.spindleInt(j, :) );
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
        function graphCut7IntensityVsTime_Surf( obj)

            times = [ 1 : length(obj.times) ] .* mean(obj.timeStep);
            distSpindle = obj.data.spindle.cut7.range;
            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_surf'); 
            surf( distSpindle, times, obj.data.spindle.cut7.spindleInt, 'LineStyle', ':', 'FaceColor', 'interp');
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
        % graphCut7MeanIntensity{{{
        function graphCut7MeanIntensity( obj)

            distSpindle = obj.data.spindle.cut7.range;
            cut7AmpAvg = mean( obj.data.spindle.cut7.spindleInt, 1);
            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_avg'); 
            plot( distSpindle, cut7AmpAvg, 'LineWidth', 2, 'Color', 'm') 
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
        % graphCut7RadialIntensity_Monopolar {{{
        function graphCut7RadialIntensity_Monopolar(obj)

            radInt = obj.data.monopolar.cut7.radialInt;
            radIntMean = mean( radInt, 1);
            radV = obj.data.monopolar.cut7.radialVal;
            f = figure('NumberTitle', 'off', 'Name', 'cut7_radial_intensity_monopolar'); 
            plot( radV, radIntMean, 'LineWidth', 2, 'Color', 'm') 
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

                    obj.graphChannelMicrotubule( )

                case 'Cut7'

                    obj.graphChannelCut7( )

            end

        end
        % }}}

        % makeMovie {{{
        function makeMovie( obj, jChannel)
            
            if find( strcmp( obj.features, 'Microtubule') ) == jChannel % if microtubule channel

                writerObj = VideoWriter([obj.path, filesep, 'mt_video.avi']);
                writerObj.FrameRate = 10;

                % open the video writer
                open(writerObj);

                % write the frames to the video
                for i=1:length(obj.data.movSpindle)

                    % convert the image to a frame
                    frame = obj.data.movSpindle(i) ;    
                    writeVideo(writerObj, frame);

                end

                % close the writer object
                close(writerObj);

            elseif find( strcmp( obj.features, 'Cut7') ) == jChannel % if microtubule channel

                writerObj = VideoWriter([obj.path, filesep, 'cut7_video.avi']);
                writerObj.FrameRate = 10;

                % open the video writer
                open(writerObj);

                % write the frames to the video
                for i=1:length(obj.data.movCut7)

                    % convert the image to a frame
                    frame = obj.data.movCut7(i) ;    
                    writeVideo(writerObj, frame);

                end

                % close the writer object
                close(writerObj);

            else

                error('unknown feature channel')

            end

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

            % Analyze if analysisData not present
            if exist( [resPathCell, filesep, 'analysisData.mat']) ~= 2

                % Initialize analysis object
                anaCell = AnalysisSingleCell( resPathCell, params.Cell.params.cellInfo.type, params.channelsToAnalyze, params.channelFeatures, params.Cell.timeStep, params.Cell.sizeVoxels);

                % analyze the cell
                anaCell.Analyze();
        
                % save this cell
                save( [resPathCell, filesep, 'analysisData.mat'], 'anaCell');

            else
                
                fprintf('      Data already exists! Loading from file...\n')
                load( [resPathCell, filesep, 'analysisData.mat'])

            end

        end
        % }}}
        
        % AnalyzeMulti {{{
        function anaCells = AnalyzeMulti( resRegExp)
            % uses AnalysisBank to analyze results from multiple movies
            % second argument can also be a cell array of multiple
            
            % Make sure we have the correct paths
            warning('off', 'MATLAB:rmpath:DirNotFound');
            rmpath( genpath(pwd) );
            warning('on', 'MATLAB:rmpath:DirNotFound');
            addpath( pwd);
            addpath( genpath( [pwd, filesep, 'classes']) );
            params = initAnalysisParams();
            addpath( genpath( params.pathParent) );

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
            
            % Compare spindle length
            %AnalysisSingleCell.GraphMulti_SpindleLength( Cells)

            % Compare Cut 7 distribution along the spindle
            %AnalysisSingleCell.GraphMulti_Cut7AlongSpindle( Cells)

            % Compare Cut 7 distribution away from pole
            AnalysisSingleCell.GraphMulti_Cut7DistributionAwayFromPole( Cells)
            
        end
        % }}}

        % GraphMulti_SpindleLength {{{
        function GraphMulti_SpindleLength( Cells)
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
                spindleLengths{jCell} = Cells{jCell}.data.spindle.mt.spindleLength;
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
        function GraphMulti_Cut7AlongSpindle( Cells)
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
                int = mean( Cells{jCell}.data.spindle.cut7.spindleInt, 1);
                if mean( int(1: floor(end/3) ) ) < mean( int( ceil(2*end/3): end) )
                    int = flip( int);
                end
                cut7int{jCell} = int/max(int);
            end
            % Distance along spindle 
            for jCell = 1 : nCells
                dist{jCell} = Cells{jCell}.data.spindle.cut7.range;
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

        % GraphMulti_Cut7DistributionAwayFromPole{{{
        function GraphMulti_Cut7DistributionAwayFromPole( Cells)
            % Compare Cut7 distibutions along the spindle

            fprintf('Comparing cut7 distirbution away from pole...\n')
            nCells = length( Cells);

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_monopolar_all'); 
            ax = axes;
            grid on; grid minor
            xlabel( 'Distance away from pole (voxels)')
            ylabel( 'Total intensity')
            set( ax, 'FontSize', 20)
            title('Cut7 intensity vs distance from pole')
            hold on;

            % Process data for plotting
            % cut7 intensity 
            for jCell = 1 : nCells
                cut7int{jCell} = Cells{jCell}.data.monopolar.cut7.radialInt;
            end

            % Distance from pole 
            for jCell = 1 : nCells
                dist{jCell} = Cells{jCell}.data.monopolar.cut7.radialVal;
            end
            
            % plot mean and shaded error
            AnalysisSingleCell.plotAreaShadedError( cut7int, dist, ax);

            for jCell = 1 : nCells
                saveas( f, [Cells{jCell}.path, filesep, 'cut7_intensity_monopolar_all.png'])
            end
            [par, ~, ~] = fileparts( Cells{1}.path );
            saveas( f, [par, filesep, 'cut7_intensity_monopolar_all.png'])
            close(f)

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
                xdata = cellfun( @(x) 1:length(x), ydata, 'UniformOutput', false)
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
                yerr(jT) = std( goodValues).^2;
            end
            max( abs( ymean - yerr) )

            % construct the extended vectors to define a region for shading
            x_vec = [x, fliplr(x)];
            y_vec = [ ymean+yerr, fliplr( ymean-yerr) ];

            % Shade the error region
            patch = fill(x_vec, y_vec, colorError);
            set(patch, 'Edgecolor', 'none');
            set(patch, 'FaceAlpha', alphaError );

            % Plot the samples
            for jSamp = 1 : nSamp 
                line( xdata{jSamp}, ydata{jSamp}, 'LineWidth', sampleWidth, 'Color', [colorSample alphaSample]) 
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
