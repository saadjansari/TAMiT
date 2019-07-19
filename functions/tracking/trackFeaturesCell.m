function trackInfo = trackFeaturesCell( featuresCell, currentCell, params);
% TrackFeatures :  Uses kalman filters to generate feature tracks given
% detected feautures for every frame. Uses the Jaqaman code for this

    if params.tracking.doTracking == 0
        disp('Tracking disabled by user via the tracking.doTracking flag in paramsInitialize.m')
        trackInfo = NaN; return
    end

    % Tracking setup {{{
    % channels available based on feature detection
    if ~isnumeric( params.config.channels) && strcmp( params.config.channels, 'all'), channelsFeatures = 1 : params.meta.numChannels;
    else, channelsFeatures = params.config.channels; end

    % channels that will be tracked
    if isnumeric( params.config.channelsToTrack) && all( ismember( params.config.channelsToTrack , channelsFeatures) ), channels = params.config.channelsToTrack;
    elseif strcmp( params.config.channelsToTrack, 'all'), channels = channelsFeatures;
    else, error('Issue with channels specification in paramsInitialize.m'), end

    % }}}

    % For each selected channel, track the detected features
    for jChannel = channels
    
        % Parse the detection data for feeding to Jaqaman tracker
        fitResults = parseDetectionData( jChannel, featuresCell(:,jChannel), params);
        
        % Load tracking parameters
        trackParams = InitModelTrackingParams( 'random_filename', params);
        
        % Call tracking function on this channel
        [ tracksFinal, kalmanInfoLink, errFlag] = trackCloseGapsKalmanSparse( ...
            fitResults, trackParams.costMatrices, trackParams.gapCloseParam, ...
            trackParams.kalmanFunctions, trackParams.probDim, ...
            trackParams.saveResults, trackParams.verbose);
        
    end    

    tracksFinal
    errFlag
    trackInfo = tracksFinal;
    
    
    % Parse Tracks
    % Initialize dynamic microtubule objects
    
    
    % Plot Tracks;



    % parseDetectionData {{{
    function fitResults = parseDetectionData( jChannel, featureInfo, params)
        
        % Temporary Code (should be addressed in the config file)
        featureType = params.config.featuresInChannels{ jChannel};
        dimTrack = featureInfo( 1).featureBank(1).dim;
%         dimTrack = 3;

        switch featureType
            case 'lines'
                for jTime = 1 : size( featureInfo, 1)

                    fitResults( jTime).probDim = 3;
                    numFeat = featureInfo(jTime).numberOfMicrotubules;
                    fitResults( jTime).num = numFeat;

                    % type of feature
                    fitResults( jTime).featureType = featureType;

                    % amplitude {{{
                    data_amp = [featureInfo(jTime).featureBank.amplitude]';
                    error_amp = zeros(numFeat, 1);
                    fitResults( jTime).amp = [ data_amp, error_amp];
                    % }}}
                    % std {{{
                    % we will use the mean std of X and Y over here

                    data_stdX = zeros(numFeat, 1);
                    data_stdY = zeros(numFeat, 1);
                    error_stdX = zeros(numFeat, 1);
                    error_stdY = zeros(numFeat, 1);
                    for jfeat = 1 : numFeat
                        data_stdX( jfeat) = featureInfo( jTime).featureBank(jfeat).fitProps.structFit.std{1}(1);
                        data_stdY( jfeat) = featureInfo( jTime).featureBank(jfeat).fitProps.structFit.std{1}(2);
                        error_stdX( jfeat) = 0.1;
                        error_stdY( jfeat) = 0.1; 
                    end
                    data_std = mean( [data_stdX, data_stdY], 2);
                    error_std = mean( [error_stdX, error_stdY], 2);

                    fitResults( jTime).std = [ data_std, error_std ];
                    % }}}
                    % start point {{{
                    % this is the MTOC
                    data_mtocX = zeros(numFeat, 1);
                    data_mtocY = zeros(numFeat, 1);
                    error_mtocX = zeros(numFeat, 1);
                    error_mtocY = zeros(numFeat, 1);
                    for jfeat = 1 : numFeat
                        data_mtocX( jfeat) = featureInfo( 1).featureBank(jfeat).mtoc(1);
                        data_mtocY( jfeat) = featureInfo( 1).featureBank(jfeat).mtoc(2);
                        error_mtocX( jfeat) = 0.1;
                        error_mtocY( jfeat) = 0.1; 
                    end
                    fitResults( jTime).xCoord = [ data_mtocX , error_mtocX ];
                    fitResults( jTime).yCoord = [ data_mtocY , error_mtocY ];
                    % }}}
                    % Poly Order {{{
                    fitResults( jTime).polyOrder = featureInfo( jTime).featureBank(1).polyOrder;
                    if dimTrack==3, fitResults( jTime).polyOrderZ = featureInfo( jTime).featureBank(1).polyOrderZ; end

                    polyCoef = zeros( numFeat, fitResults(jTime).polyOrder+1, 3);
                    
                    for jmt = 1 : numFeat
                        polyCoef( jmt, :, 1) = featureInfo( jTime).featureBank(jmt).fitCoef{1};
                        polyCoef( jmt, :, 2) = featureInfo( jTime).featureBank(jmt).fitCoef{2};
                        if dimTrack==3
                            coefZ = nan( 1, fitResults(jTime).polyOrder);
                            startIdx = fitResults(jTime).polyOrder - fitResults(jTime).polyOrderZ + 1;
                            polyCoef( jmt, startIdx:end, 3) = featureInfo( jTime).featureBank(jmt).fitCoef{3};
                        end
                    end
                    fitResults( jTime).polyCoef = polyCoef;
                    % }}}

                    disableThis = 1;
                    if ~disableThis
                    % length {{{
                    % find length of each microtubule
                    lengthmt = zeros( numFeat, 1 );
                    for jmt = 1 : length( featureInfo(jTime).featureBank)
                        lengthmt(jmt) = featureInfo(jTime).featureBank(jmt).findLength;
                    end
                    % Find Way to evaluate the error in lengths.
                    fitResults( jTime).length = [ lengthmt , zeros(numFeat, 1)];

                    % }}}
                    % theta {{{
                    % find mean orientation of each microtubule
                    thetamt = zeros( numFeat, 1);
                    for jmt = 1 : length( featureInfo(jTime).featureBank)
                        [thetamt(jmt), ~] = featureInfo(jTime).featureBank(jmt).findOrientationMean;
                    end
                    fitResults( jTime).theta = [ thetamt, zeros(numFeat, 1) ];

                    % Find a way to evaluate the error in the mean orientation

                    % }}}
                    end

                end
            case 'spots'
                for jTime = 1 : size( featureInfo, 1)

                    fitResults( jTime).probDim = 3;
                    numFeat = length( featureInfo(jTime).numberOfSpots);
                    fitResults( jTime).num = numFeat;

                    % type of feature
                    fitResults( jTime).featureType = featureType;

                    % amplitude {{{
                    fitResults( jTime).amp = [ [featureInfo(jTime).featureBank.amplitude]' , zeros(numFeat, 1) ];
                    % }}}
                    % std {{{
                    % we will use the mean std of X and Y over here
                    stdX = [featureInfo(jTime).featureBank.std{1}(1)]';
                    stdY = [featureInfo(jTime).featureBank.std{1}(2)]';
                    fitResults( jTime).std = [ mean( [stdX; stdY], 1)', zeros(numFeat, 1) ];
                    % }}}
                    % start point {{{
                    % this is the MTOC
                    mtocX = [featureInfo(jTime).featureBank.coord(1)]';
                    mtocY = [featureInfo(jTime).featureBank.coord(2)]';
                    try; mtocZ = [featureInfo(jTime).featureBank.coord(3)]'; end
                    fitResults( jTime).xCoord = [mtocX , zeros(numFeat, 1) ];
                    fitResults( jTime).yCoord = [mtocY , zeros(numFeat, 1) ];
                    try; fitResults( jTime).zCoord = [mtocZ , zeros(numFeat, 1) ]; end
                    % }}}
                    % Poly Order {{{
                    fitResults( jTime).polyOrder = 0;
                    fitResults( jTime).polyCoef = NaN;
                    % }}}
                end
                error('Spots tracking not set up as of yet')
            otherwise
                error('What is this mystery feature my master? I have been slayed.')
        end

    end
    % }}}

end

