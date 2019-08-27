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
        mov
        data
        folderName
    end

    methods (Access = public )

        % AnalsisSingleCell {{{
        function obj = AnalysisSingleCell( pathResCell, cellType, channels, features, timeStep, sizeVoxels, mov)
            % Contruct the Analysis Object
            
            obj.path = pathResCell;
            if ~strcmp(cellType, 'Mitosis') && ~strcmp( cellType, 'Interphase'),
                error('AnalysisBank: unknown cellType, allowed values are Mitosis and Interphase.')
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

            if nargin > 6
                obj.mov = mov;
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
                fprintf('   Analyzing channel %d with feature %s\n', channel, feature); 

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
            % FIXME cut7 does not have data, so we'll  use microtubule data
            if strcmp( obj.features{jChannel}, 'Cut7')
                % find MT channel
                refChannel = obj.channels( find( strcmp( obj.features, 'Microtubule') ) );
                obj.times = obj.findTimes( refChannel);
            end 

            for jFrame = 1 : length( obj.times)

                % Current time frame
                jTime = obj.times( jFrame);

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

            % construct data file name and load it
            timeDataFile = ['C' num2str( channel) '_T' num2str( time) '_global.mat']; 

            % Load data file
            
            % FIXME cut7 does not have data, so we'll load up microtubule data for it to get spindle lengths
            if strcmp( obj.features{jChannel}, 'Cut7')
                % find MT channel
                refChannel = obj.channels( find( strcmp( obj.features, 'Microtubule') ) );
                timeDataFile = ['C' num2str( refChannel) '_T' num2str( time) '_global.mat']; 
            end 
            timeData = load( timeDataFile);

            % Initialize main feature from loaded data
            eval( ['mainFeature = ' timeData.featureMainStruct.type '.loadFromStruct( timeData.featureMainStruct);'] );
            %if strcmp( obj.features{jChannel}, 'Cut7')
                %timeData.featureMainStruct.type
                %mainFeature
            %end 

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
                    obj.data.spindleLength( jTime) = obj.analyzeSpindle( mainFeature);
                case 'Monopolar'
                    error('analyzeFeatureMicrotubule: Monopolar still in development')
                case 'InterphaseBank'
                    error('analyzeFeatureMicrotubule: InterphaseBank still in development')
                otherwise
                    error('analyzeFeatureMicrotubule: unknown mainFeature type')
            end

            % Capture frame for movie
            img = obj.mov( :,:,:, obj.times(jTime), 2);
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
        function spindleLength = analyzeSpindle( obj, mainFeature)
            % analyze a microtubule spindle            

            % Find spindle lengths
            spindleLength = norm( (mainFeature.featureList{1}.startPosition - mainFeature.featureList{1}.endPosition) .* obj.sizeVoxels );

            % Any other things?

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

            end

        end
        % }}}
        % graphSpindleLength {{{
        function graphSpindleLength( obj)
            % graph the spindle length vs time

            times = [ 1 : length(obj.times) ]* mean(obj.timeStep);
            % plot spindle Length
            f = figure('NumberTitle', 'off', 'Name', 'spindle_length'); 
            plot( times, obj.data.spindleLength, 'LineWidth', 2, 'Color', 'm') 
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
       % }}} 
       
        % Cut7 {{{
        % analyzeFeatureCut7 {{{
        function analyzeFeatureCut7( obj, mainFeature, jTime)
            % do cut7 analysis
            
            cTime = obj.times( jTime);

            % get spindle start, spindle end from microtubule channel
            startPos = mainFeature.featureList{1}.startPosition;
            endPos = mainFeature.featureList{1}.endPosition;

            % get cut7 frame 
            channelCut7 = obj.channels( find( strcmp( obj.features, 'Cut7') ) );
            cut7Frame = obj.mov(:,:,:, cTime, channelCut7);
            % take a max z projection
            cut7Frame = squeeze( max( cut7Frame, [], 3) );
            cut7amp = Cell.findAmplitudeAlongLine( cut7Frame, startPos(1:2), endPos(1:2)); 

            % normalize it from 0 to 1
            normRange = 0:0.02:1;
            if ~isfield( obj.data, 'cut7Range');
                obj.data.cut7Range = normRange;
            end

            % Store amplitude data normalized against spindle length
            obj.data.cut7Amp(jTime, :)=interp1( linspace(0,1,length( cut7amp) ), cut7amp, normRange);

        end
        % }}}
        % graphChannelCut7 {{{
        function graphChannelCut7( obj)

            % graph cut7 intensity along the spindle averaged over time
            obj.graphCut7MeanIntensity();
            
            % graph cut7 intensity along the spindle with 3rd dimension corresponding to time 
            obj.graphCut7IntensityVsTime_Overlay();
            obj.graphCut7IntensityVsTime_Surf();

        end
        % }}}
        % graphCut7IntensityVsTime_Overlay {{{
        function graphCut7IntensityVsTime_Overlay( obj)

            times = [ 1 : length(obj.times) ] .* mean(obj.timeStep);
            distSpindle = obj.data.cut7Range;

            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_overlay'); 
            hold on;
            for j = 1 : length(times)
                plot( distSpindle, obj.data.cut7Amp(j, :) );
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
            distSpindle = obj.data.cut7Range;
            f = figure('NumberTitle', 'off', 'Name', 'cut7_intensity_spindle_time_surf'); 
            surf( distSpindle, times, obj.data.cut7Amp, 'LineStyle', ':', 'FaceColor', 'interp');
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

            distSpindle = obj.data.cut7Range;
            cut7AmpAvg = mean( obj.data.cut7Amp, 1);
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
            
            % FIXME
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

            end

        end
        % }}}
        
    end

    methods ( Static = true, Access = public )

        % AnalyzeSingle {{{
        function dataCell = AnalyzeSingle( resPathCell, params)
            % uses AnalysisBank to analyze results from a single movie
            
            % Make sure we have the correct paths
            warning('off', 'MATLAB:rmpath:DirNotFound');
            rmpath( genpath(pwd) );
            warning('on', 'MATLAB:rmpath:DirNotFound');
            addpath( pwd);
            addpath( [pwd, filesep, 'functions', filesep, 'analysis'] );
            params = initAnalysisParams();
            addpath( genpath( params.pathParent) );
            addpath( [pwd, filesep, 'classes/']);
            addpath( genpath( [pwd, filesep, 'functions/external']) );

            % if params not provided, find the params
            if nargin < 2
                params = initAnalysisParams();
            end

            % check if resPathCell is a fullpath, if not then make it full
            [pth, ~, ~] = fileparts( resPathCell);
            if isempty( pth)
                resPathCell = fullfile( params.pathParent, resPathCell);
            end
            fprintf( 'resultsFolderPath = %s\n', resPathCell);
           
            params.Cell = load( [resPathCell, filesep, 'params.mat']);
            params.resPath = resPathCell;

            % TEMPORARY. IMPORTING THE MOVIE AND LOADING SOME INFORMATION.
            % THIS SHOULD BE DONE BY SINGLECELL
            movInfo = AnalysisSingleCell.importMovie( params);

            % Initialize analysis object
            analysisObj = AnalysisSingleCell( resPathCell, params.cellType, params.channelsToAnalyze, params.channelFeatures, movInfo.timeStep, movInfo.sizeVoxels, movInfo.mov);

            % analyze the cell
            analysisObj.Analyze();
        
        end
        % }}}
        
        % AnalyzeMulti {{{
        function AnalyzeMulti( paramsPath, resRegExp)
            % uses AnalysisBank to analyze results from multiple movies
            % second argument can also be a cell array of multiple
            
            % Make sure we have the correct paths
            warning('off', 'MATLAB:rmpath:DirNotFound');
            rmpath( genpath(pwd) );
            warning('on', 'MATLAB:rmpath:DirNotFound');
            addpath( pwd);
            addpath( [pwd, filesep, 'functions', filesep, 'analysis'] );
            params = initAnalysisParams();
            addpath( genpath( params.pathParent) );
%             addpath( [pwd, filesep, 'classes/']);

            % Get all folder names in parent
            f = dir( params.pathParent);
            folds = folds.name;
            folds = folds( f.isdir);

            % Remove all folders starting with dot 
            folds = folds( cellfun( @isempty, regexp( folds, '^[.]') ) );

            % Ensure that we only have the results folders in the correct format (start with YYMMDD_HHMM)
            folds = folds( cellfun( @(x) ~isempty(x), regexp( folds, '^\d{6}_\d{4}') ) );
            
            % match resRegExp
            tempExp = '1095_50msG_100msR_7Z';
            fprintf('Note in AnalysisBank.AnalyzeMulti: temporary regexpression defined in code')
            folds = folds( cellfun( @(x) ~isempty(x), regexp( folds, tempExp) ) );
            numCells = length( folds);
            fprintf('AnalysisBank.AnalyzeMulti found %d cells to analysis', numCells);

            % for each cell we'll run an AnalyzeSingle function
            for jCell = 1 : numCells
                dataCell = AnalysisBank.AnalyzeSingle( folds( jCell), params ); 
            end

        end
        % }}}

        % importMovie {{{
        function dataCell = importMovie( params)
            % imports the movie and loads up important data about the experiment

            if params.importedFromSummit 
                moviePathSummit = params.Cell.params.cellinfo.moviePath;
                moviePath = erase( moviePathSummit, params.rmStringSummit);
                moviePath = [params.addString, moviePath];
            else
                moviePath = params.Cell.params.cellInfo.moviePath;
            end

            % import the movie
            addpath( 'functions/')
            cellData = importSingleCell( moviePath); 

            % load the cell movie, voxel sizes, and calculate time steps
            dataCell.mov = cellData.cell3D;
            dataCell.sizeVoxels = [ cellData.metaData.sizeVoxelsX, cellData.metaData.sizeVoxelsY, cellData.metaData.sizeVoxelsZ];

            % find mean time for all z-slices
            timesTC = squeeze(mean( cellData.planeTimes, 1) ); 
            timeSteps = diff(timesTC);

            % remove any nans and then use the median timeStep as the timeStep
            timeSteps( isnan( timeSteps) ) = 0;
            dataCell.timeStep = median( timeSteps);
%             fprintf('Time Step = %4.2f sec\n', dataCell.timeStep)

        end
        % }}}

    end

end
