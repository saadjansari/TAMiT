function trackInfo = trackFeaturesCell( mov, featuresCell, currentCell, params);
% TrackFeatures :  Uses kalman filters to generate feature tracks given
% detected feautures for every frame. Uses the Jaqaman code for this

    if params.tracking.doTracking == 0
        disp('Tracking disabled by user via the tracking.doTracking flag in paramsInitialize.m')
        trackInfo = NaN; return
    end

    % Tracking setup {{{

    % Specify for which channels to track the detected features.
    % channels available based on feature detection
    if ~isnumeric( params.channels) && strcmp( params.channels, 'all'), channelsFeatures = 1 : params.meta.numChannels;
    else, channelsFeatures = params.channels; end

    % channels that will be tracking
    if isnumeric( params.tracking.channels) && all( ismember( params.tracking.channels , channelsFeatures) ), channels = params.tracking.channels;
    elseif strcmp( params.tracking.channels, 'all'), channels = channelsFeatures;
    else, error('Issue with channels specification in paramsInitialize.m'), end

    % }}}

    % For each selected channel, track the detected features
    for jChannel = trackChannels
    
        % Parse the detection data for feeding to Jaqaman tracker
        fitResults = parseDetectionData( jChannel, featuresCell(:,jChannel), params);
        
        % Call tracking function on this channel
        [ tracksFinal, kalmanInfoLink, errFlag] = trackCloseGapsKalmanSparse( ...
            fitResults, trackParams.costMatrices, trackParams.gapCloseParam, ...
            trackParams.kalmanFunctions, trackParams.probDim, ...
            trackParams.saveResults, trackParams.verbose);
        
        % Save tracking results
        trackResults.Channel( jChannel).tracksFinal = tracksFinal;
        trackResults.Channel( jChannel).kalmanInfoLink = kalmanInfoLink;
        trackResults.Channel( jChannel).errFlag = errFlag;
        
        % Overlay tracked features
        mov = plotTracks( Movies.Channel{ jChannel }, movieInfo.fitResults( 1).frame, tracksFinal );
        
        h = figure;
        movie( h, mov, 1, 0.5)
%         % Overlay tracked features
%         savedir = [pwd '\'];
%         filenamebase = Movies.MetaData.fileName;

%         % Save MIP images
%         firstImageFile = [savedir,'Results\test_1.tif'];
%         for jFrame = 1 : length(fitResults)
%             imPlane = mat2gray( max( Movies.Cell( 1).Channel{2 } (:, :, :, jFrame), [], 3) );
%             if jFrame == 1
%                 imwrite( imPlane, firstImageFile)
%             else
%                 ImageFileName = [savedir,'Results\test_', num2str(jFrame), '.tif'];
%                 imwrite( imPlane, ImageFileName, 'WriteMode', 'append')
%             end
%         end
%         
%         movieName = [filenamebase,'_movie']; %'testND(wt)_fit_movie'; 
% %         firstImageFile = [savedir,filenamebase,'_cell_' num2str(1) '_fit_frame_' num2str(1) '.png']; % 'testND(wt)_fit_1.png';
%         minTrackLength = 3; 
%         dragtail_length = 0;
%         Imin = 0;    %min(Movies.Cell(current_cell).Channel{1}(:));
%         Imax = max(Movies.Cell( 1).Channel{2 }(:));
% 
%         % Overlay lines and determine lengths
%         overlayTracksMovieNewMB(tracksFinal,[1 190],dragtail_length,...
%             1,movieName,0,[],[],2,[],[],[],[],3,1,firstImageFile,savedir,minTrackLength,...
%             [],'mov',[Imin 0.55*Imax]);
%         
%         overlayTracksMovieNewMB(tracksFinal,startend,dragtailLength,...
%             saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
%             imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
%             colorTracks,firstImageFile,dir2saveMovie,minLength,plotFullScreen,...
%             movieType,intensityRange)
        
%         overlayTracksMovieKC3D(tracksFinal,movieInfo,startend,dragtailLength,...
%             saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
%             imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
%             colorTracks,firstImageFile,dir2saveMovie,minLength,plotFullScreen,...
%             movieType,intensityRange,original_image,curr_cell)
        
    end    

    % parseDetectionData {{{
    function fitResults = parseDetectionData( jChannel, featureInfo, params)
        
        % Temporary Code (should be addressed in the config file)
        if jChannel==1, featureType='lines';
        elseif jChannel==2, featureType='spots';
        end

        switch featureType
            case 'lines'
                for jTime = 1 : size( featureInfo, 1)
                    % amplitude {{{
                    fitResults( jTime).amp = [featureInfo(jTime).featureBank.amplitude]
                    % }}}
                    % std {{{
                    % we will use the mean std of X and Y over here
                    stdX = [featureInfo(jTime).featureBank.std{1}(1)]
                    stdY = [featureInfo(jTime).featureBank.std{1}(2)]
                    fitResults( jTime).std = mean( [stdX, stdY] );
                    % }}}
                    % start point {{{
                    % this is the MTOC
                    mtocX = [featureInfo(jTime).featureBank.mtoc(1)]
                    mtocY = [featureInfo(jTime).featureBank.mtoc(2)]
                    mtocZ = [featureInfo(jTime).featureBank.mtoc(3)]
                    fitResults( jTime).mtoc = [mtocX ; mtocY ; mtocZ ];
                    % }}}
                    % length {{{
                    % find length of each microtubule
                    lengthmt = zeros( 1, length( featureInfo(jTime).featureBank) )
                    for jmt = 1 : length( featureInfo(jTime).featureBank)
                        lengthmt(jmt) = featureInfo(jTime).featureBank(jmt).findLength;
                    end
                    fitResults( jTime).length = lengthmt;
                    % }}}
                    % theta {{{
                    % find mean orientation of each microtubule
                    thetamt = zeros( 1, length( featureInfo(jTime).featureBank) )
                    for jmt = 1 : length( featureInfo(jTime).featureBank)
                        [thetamt(jmt), ~] = featureInfo(jTime).featureBank(jmt).findOrientationMean;
                    end
                    fitResults( jTime).theta = thetamt;
                    % }}}
                end
            case 'spots'

            otherwise
                error('What is this mystery feature oh master? I have been slayed.')
        end

    end
    % }}}
    
end

