function [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
    nnDistFeatures,prevCost,errFlag] = linkFeaturesKalmanSparseSA(movieInfo,...
    costMatName,costMatParam,kalmanFunctions,probDim,filterInfoPrev,...
    prevCost,verbose)
%LINKFEATURESKALMAN links features between consecutive frames using LAP and possibly motion propagation using the Kalman filter
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
%    nnDistFeatures,prevCost,errFlag] = linkFeaturesKalmanSparse(movieInfo,...
%    costMatName,costMatParam,kalmanFunctions,probDim,filterInfoPrev,...
%    prevCost,verbose)
%
%INPUT  movieInfo      : Array of size equal to the number of frames
%                        in a movie, containing the fields:
%             .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .zCoord      : z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if problem is 2D. Default: zeros.
%             .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%             .num         : Number of features. 
%                            Optional. Calculated if not supplied.  
%             .allCoord    : x,dx,y,dy,[z,dz] of features collected in one
%                            matrix. Optional. Calculated if not supplied.
%             .nnDist      : Distance from each feature to its nearest
%                            neighbor. Optional. Calculated if not supplied.
%       costMatName    : Name of cost matrix function used for linking.
%       costMatParam   : Parameters needed for cost matrix calculation. 
%                        Structure with fields specified by particular
%                        cost matrix used (costMatName).
%       kalmanFunctions: Names of Kalman filter functions for self-adaptive
%                        tracking. Structure with fields:
%             .reserveMem   : Reserves memory for kalmanFilterInfo.
%             .initialize   : Initializes the Kalman filter for an appearing
%                             feature.
%             .calcGain     : Calculates the Kalman gain after linking.
%                        For non-self-adaptive tracking, enter [].
%                        Optional. Default: [].
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                        Optional. If not input, dimensionality will be
%                        derived from movieInfo.
%       filterInfoPrev : Structure array with number of entries equal to
%                        number of frames in movie. Contains at least the
%                        fields:
%             .stateVec    : Kalman filter state vector for each feature in frame.
%             .stateCov    : Kalman filter state covariance matrix for each feature in frame.
%             .noiseVar    : Variance of state noise for each feature in frame.
%                            Optional. Enter [] or nothing if not to be used.
%       prevCost       : Matrix of costs of actual assignments in previous
%                        round of linking. Optional. Default: Empty.
%       verbose        : 1 to show calculation progress, 0 otherwise.
%                        Optional. Default: 1.
%
%OUTPUT trackedFeatureIndx: Connectivity matrix of features between frames.
%                           Rows indicate continuous tracks, while columns 
%                           indicate frames. A track that ends before the
%                           last frame is followed by zeros, and a track
%                           that starts in a frame after the first frame
%                           is preceded by zeros. 
%       trackedFeatureInfo: The positions and "intensities" of the tracked
%                           features, in the same units as input. 
%                           Number of rows = number of tracks.
%                           Number of columns = 8*number of frames.
%                           Each row consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           NaN is used to indicate frames where the tracks
%                           do not exist.
%       kalmanFilterInfo  : Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields
%                           defined in kalmanFunctions.reserveMem (at
%                           least stateVec, stateCov and noiseVar).
%       nnDistFeatures    : Matrix indicating the nearest neighbor
%                           distances of features linked together within
%                           tracks.
%       prevCost          : Matrix of costs of actual assignments.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must not be empty.
%
%Khuloud Jaqaman, March 2007

%% Output

trackedFeatureIndx = [];
trackedFeatureInfo = [];
kalmanFilterInfo = [];
nnDistFeatures = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--linkFeaturesKalmanSparse: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether tracking is self-adaptive 
if nargin < 4 || isempty(kalmanFunctions)
    kalmanFunctions = [];
    selfAdaptive = 0;
else
    selfAdaptive = 1;
end

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else    
    probDimT = 2;
end

%assign problem dimensionality if not input
if nargin < 5 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--linkFeaturesKalmanSparse: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%check whether a priori Kalman filter information is given
if nargin < 6 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%check whether previous costs have been input
if nargin < 7 || isempty(prevCost)
    prevCost = [];
end

%check whether verbose
if nargin < 8 || isempty(verbose)
    verbose = 1;
end

%exit if there are problems with input
if errFlag
    disp('--linkFeaturesKalmanSparse: Please fix input parameters.');
    return
end

%% preamble

%get number of frames in movie
numFrames = length(movieInfo);

%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iFrame = 1 : numFrames
        movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    if probDim == 2
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord];
        end
    else
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord movieInfo(iFrame).zCoord];
        end
    end
end

%calculate nearest neighbor distance for each feature in each frame
if ~isfield(movieInfo,'nnDist')

    for iFrame = 1 : numFrames
        
        switch movieInfo(iFrame).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there is more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iFrame).allCoord(:,1:2:end),...
                    movieInfo(iFrame).allCoord(:,1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);

        end

        %store nearest neighbor distance
        movieInfo(iFrame).nnDist = nnDist;

    end

end

%% Linking

% get number of features in first frame and second frame
NParticles = 6;
N1 = movieInfo(1).num;
N2 = movieInfo(2).num;

%make an array of the number of features per frame
numFeatures = zeros(numFrames,1);
for iFrame = 1 : numFrames
    if iFrame == 1
        numFeatures( iFrame) = max( [ N1 , N2 ] );
    else
        numFeatures(iFrame) = NParticles;
    end
end

%reserve memory for kalmanFilterInfo
if selfAdaptive

    % -- USER DEFINED FUNCTION -- %
    eval(['kalmanFilterInfo = ' kalmanFunctions.reserveMem ...
        '(numFrames,numFeatures,probDim);']);
    
else
    
    kalmanFilterInfo = zeros(numFrames,1);

end

%fill the feature indices in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%fill the nearest neighbor distances of features in first frame
nnDistFeatures = movieInfo(1).nnDist;

% %initialize Kalman filter for features in 1st frame
% if selfAdaptive
% 
%     % -- USER DEFINED FUNCTION -- %
%     if usePriorInfo %use a priori information if available
%         kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
%         kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
%         kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
%     else
%         eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
%             '(movieInfo(1),probDim,costMatParam);'])
%         kalmanFilterInfo(1).stateVec = filterInit.stateVec;
%         kalmanFilterInfo(1).stateCov = filterInit.stateCov;
%         kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
%     end
%     
% end

%store the costs of previous links for features in first frame
%if no previous costs have been input
%in this case, store NaN since there are no previous links
if isempty(prevCost)
    prevCost = NaN(movieInfo(1).num,1);
else
    prevCost = max(prevCost(:))*ones(movieInfo(1).num,1);
end
prevCostStruct.all = prevCost;
prevCostStruct.max = max(prevCost(:));
prevCostStruct.allAux = [];

%assign the lifetime of features in first frame
featLifetime = ones(movieInfo(1).num,1);

% % % %for paper - get number of potential link per feature
% % % numPotLinksPerFeature = [];

%get number of particles in whole movie and calculate a worst-case scenario
%number of tracks
%it can be that the final number of tracks is even larger than this worst
%case scenario. Every time the auxiliary matrices (defined below) run out
%of rows, another "numTracksWorstCase" rows are added to them.
if round(sum(numFeatures)/10) == 0
    numTracksWorstCase = 1;
else
    numTracksWorstCase = round(sum(numFeatures)/10);
end

%initialize auxiliary matrices for storing information related to tracks
%that end in the middle of the movie
trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
nnDistFeaturesAux = NaN(numTracksWorstCase,numFrames);
prevCostAux = NaN(numTracksWorstCase,numFrames);
rowEnd = numTracksWorstCase;

%initialize progress display
if verbose
    progressText(0,'Linking frame-to-frame');
end

% -------- SA edit --------

%initialize Kalman filter for features in 1st frame
% if selfAdaptive
% 
%     % -- USER DEFINED FUNCTION -- %
%     if usePriorInfo %use a priori information if available
%         kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
%         kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
%         kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
%     else
%         eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
%             '(movieInfo(1),probDim,costMatParam);'])
%         kalmanFilterInfo(1).stateVec = filterInit.stateVec;
%         kalmanFilterInfo(1).stateCov = filterInit.stateCov;
%         kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
%     end
%     
% end

% Main Premise : There are always six particles in our movie and they
% usually contribute to the detected amplitude of our spots. They may not
% be identified as unique spots, but their intensity is picked up and
% included in the amplitude of one that is picked up.

% In this part, we will use the terminology of shadow particles and
% concealing particles. Say we detect 4 unique spots in frame 1, and 6
% in frame 2. We will say that there are 6 spots in frame 1, and we can
% only see 4 because 2 are being concealed by the 4 spots that we see. For
% each shadow particle, we'll assign an identity to it that is a copy of
% the particle that is concealing it. To determine the concealing particle,
% we'll pick one that minimizes the cost of linking (distance and amplitude
% considerations)

% go over all frames
for iFrame = 1 : numFrames-1

    %get number of features
    numFeaturesFrame1 = movieInfo(iFrame).num; %in 1st frame
    numFeaturesFrame2 = movieInfo(iFrame+1).num; %in 2nd frame
    
    % first job is to match the number of particles between frame 1 and 2 
    if iFrame == 1
       
        
        % find max number of particles, number of shadow particles and the
        % concealing frame
        [ numTrue, ~] = max( [ numFeaturesFrame1, numFeaturesFrame2 ] );
        numShadows = abs( numFeaturesFrame1 - numFeaturesFrame2 );
        [ ~, suspectFrame] = min( [ numFeaturesFrame1, numFeaturesFrame2 ] );
        suspectFrame = iFrame + suspectFrame - 1;
        if suspectFrame == iFrame ; refFrame = iFrame + 1; else refFrame = iFrame; end
        
        numConcealers = movieInfo( suspectFrame).num;
        
        % We will run a loop that assigns a concealer to each shadow
        % particle. We will test out all combinations, and pick the one
        % that minimizes the linking cost.
        % initialize possible combinations
        
        permutes = permn( 1 : numConcealers, numShadows);
        % each column is a shadow, each row is a possible combination, and
        % each each value is the index of the concealing particle.
        matches = cell( 1, size( permutes, 1) );
        costEff = zeros( 1, size( permutes, 1) );
        suspectInfoCell = cell(1, size( permutes, 1) );
        for jComb = 1 : size( permutes, 1)
            
            % extract the relevant combination
            comb = permutes( jComb, :);
            
            suspectInfo = movieInfo( suspectFrame);
            refInfo = movieInfo( refFrame);
            
            % divide amplitude intensity appropriately
            for jSpot = unique( comb)
                 numConceals = 1 + sum( comb == jSpot);
                 suspectInfo.SpotAmp( jSpot) = suspectInfo.SpotAmp( jSpot) / numConceals;
            end
                
            % create an identity for the victims
            for jVictim = 1 : numShadows
                
                suspectInfo.allCoord = [ suspectInfo.allCoord ; suspectInfo.allCoord( comb( jVictim), :) ];
                suspectInfo.SpotAmp = [ suspectInfo.SpotAmp, suspectInfo.SpotAmp( comb( jVictim) ) ];
                suspectInfo.SpotAmpErr = [ suspectInfo.SpotAmpErr, suspectInfo.SpotAmpErr( comb( jVictim) ) ];
                suspectInfo.SpotSigmaXY = [ suspectInfo.SpotSigmaXY, suspectInfo.SpotSigmaXY( comb( jVictim) ) ];
                suspectInfo.SpotSigmaZ = [ suspectInfo.SpotSigmaZ, suspectInfo.SpotSigmaZ( comb( jVictim) ) ];
                suspectInfo.SpotSigmaXYErr = [ suspectInfo.SpotSigmaXYErr, suspectInfo.SpotSigmaXYErr( comb( jVictim) ) ];
                suspectInfo.SpotSigmaZErr = [ suspectInfo.SpotSigmaZErr, suspectInfo.SpotSigmaZErr( comb( jVictim) ) ];
                suspectInfo.num = numTrue;
                
            end
            
            suspectInfoCell{ jComb} = suspectInfo;
            
            % calculate amplitude bias for cost matrix calculation
            biasAmp = zeros( 1, numShadows);
            for jVictim = 1 : numShadows
                
                ampC = suspectInfo.SpotAmp( comb( jVictim) );
                ampO = suspectInfo.SpotAmp; ampO ( comb( jVictim) ) = []; ampO = mean( ampO);
                
                biasAmp( jVictim) = 1 + 0.5 * abs( ampC - ampO ) / ampO; 
                
            end
            
            % find cost matrix from nearest neighbor distance calculation
            nnDist = createDistanceMatrix( suspectInfo.allCoord(:,1:2:end),...
            refInfo.allCoord(:,1:2:end));

            % add bias to nnDist
            for jVictim = 1 : numShadows
            
                nnDist( numConcealers + jVictim, :) = biasAmp( jVictim) .* nnDist( numConcealers + jVictim, :);
                
            end

            % find the effective minimum linking cost for unique matchings
            % Will have to try all combinations of unique matches, for best
            % global match
            pp = perms( 1 : numTrue);
            costs = zeros( size( pp,1), 1);
            for jMatch = 1 : size( pp, 1)
                
                idx2 = pp( jMatch, :);
                idx1 = 1 : length( idx2);
                
                costs( jMatch) = sum( nnDist( sub2ind(size(nnDist), idx1, idx2) ) );
                
            end
            
            % pick the combination with the min overall cost
            [costEff( jComb), idx] = min( costs);
            matches{ jComb} = [ idx1' , [ pp( idx, :)]' ];
                
        end
        
        % We have our first linking assignment between frame1 and frame 2
        [~, idx] = min( costEff);
        match = matches{ idx};
        
        link12 = match(:, 2);
        [~, idxx] = sort( link12, 'ascend');
        link21 = match( idxx,1);
        
        % Store identity of victim in movieInfo2
        movieInfo2( suspectFrame) = suspectInfoCell{ idx};
        movieInfo2( refFrame) = movieInfo( refFrame);

        % Now we will increase the number of particles in frame 2 to be
        % equal to what we expect.
        numFrame2 = NParticles;
        
        numShadows = abs( numFrame2 - numTrue );
        initFrame = iFrame;
        suspectFrame = iFrame + 1;
        
        numConcealers = movieInfo2( suspectFrame).num;
        numVictims = NParticles - numConcealers;
        
        % use amplitude consideration to convert numConcealers to 6
        % particles
        ampSpot = sum( movieInfo2(suspectFrame).SpotAmp) / NParticles;
        % store number of duplicates
        concealed = zeros(1, numConcealers);
        
        for jCon = 1 : numConcealers
           
            concealed( jCon) = movieInfo2(suspectFrame).SpotAmp( jCon) / ampSpot;
            
        end
        
        concealedDef = floor( concealed);
        [ temp, ind] = sort( mod( concealed, 1), 'descend');
        concealedDef( ind( 1 : NParticles - sum(concealedDef) ) ) = ...
            concealedDef( ind( 1 : NParticles - sum(concealedDef) ) ) + 1;
        
        concealedDef = concealedDef -1;
        
        % Assign identities to the victim particles.
        for jCon = 1 : numConcealers
            
            if concealedDef( jCon) >= 1
               
                nV = concealedDef(jCon);
                amp = movieInfo2( suspectFrame).SpotAmp( jCon) / (nV+1) ;
                movieInfo2( suspectFrame).SpotAmp( jCon) = amp;
                movieInfo2( suspectFrame).SpotAmp = movieInfo2( suspectFrame).SpotAmp( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).SpotAmpErr = movieInfo2( suspectFrame).SpotAmpErr( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).SpotSigmaXY = movieInfo2( suspectFrame).SpotSigmaXY( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).SpotSigmaZ = movieInfo2( suspectFrame).SpotSigmaZ( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).SpotSigmaXYErr = movieInfo2( suspectFrame).SpotSigmaXYErr( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).SpotSigmaZErr = movieInfo2( suspectFrame).SpotSigmaZErr( [1:end, jCon*ones( 1, nV) ] );
                movieInfo2( suspectFrame).allCoord = [ movieInfo2( suspectFrame).allCoord ; repmat( movieInfo2( suspectFrame).allCoord( jCon, :), nV, 1) ];
                movieInfo2( suspectFrame).num = NParticles;

            end
            
        end
        
                %initialize Kalman filter for features in 1st frame
        if selfAdaptive

            % -- USER DEFINED FUNCTION -- %
            if usePriorInfo %use a priori information if available
                kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
                kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
                kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
            else
                eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
                    '(movieInfo2(1),probDim,costMatParam);'])
                kalmanFilterInfo(1).stateVec = filterInit.stateVec;
                kalmanFilterInfo(1).stateCov = filterInit.stateCov;
                kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
            end

        end
        
    else
        % not the first frame. Here we assume that we are starting with a
        % fixed number of particles, say 6. The second frame may or may not
        % have all 6 of these particles (but must be less than 6)
        
        % Frame i has NParticles. Frame i+1 may or mat not have N
        % particles.
        % Extract number of particles in frame i+1
        numConcealers = numFeaturesFrame2;
        numVictims = NParticles - numConcealers;
        suspectFrame = iFrame+1;
        refFrame = iFrame;
        
        % We will run a loop that assigns a concealer to each shadow
        % particle. We will test out all combinations, and pick the one
        % that minimizes the linking cost.
        % initialize possible combinations
        
        permutes = permn( 1 : numConcealers, numVictims);
        % each column is a shadow, each row is a possible combination, and
        % each each value is the index of the concealing particle.
        matches = cell( 1, size( permutes, 1) );
        costEff = zeros( 1, size( permutes, 1) );
        suspectInfoCell = cell(1, size( permutes, 1) );
        for jComb = 1 : size( permutes, 1)
            
            % extract the relevant combination
            comb = permutes( jComb, :);
            
            suspectInfo = movieInfo( suspectFrame);
            refInfo = movieInfo2( refFrame);
            
            % divide amplitude intensity appropriately
            for jSpot = unique( comb)
                 numConceals = 1 + sum( comb == jSpot);
                 suspectInfo.SpotAmp( jSpot) = suspectInfo.SpotAmp( jSpot) / numConceals;
            end
                
            % create an identity for the victims
            for jVictim = 1 : numShadows
                
                suspectInfo.allCoord = [ suspectInfo.allCoord ; suspectInfo.allCoord( comb( jVictim), :) ];
                suspectInfo.SpotAmp = [ suspectInfo.SpotAmp, suspectInfo.SpotAmp( comb( jVictim) ) ];
                suspectInfo.SpotAmpErr = [ suspectInfo.SpotAmpErr, suspectInfo.SpotAmpErr( comb( jVictim) ) ];
                suspectInfo.SpotSigmaXY = [ suspectInfo.SpotSigmaXY, suspectInfo.SpotSigmaXY( comb( jVictim) ) ];
                suspectInfo.SpotSigmaZ = [ suspectInfo.SpotSigmaZ, suspectInfo.SpotSigmaZ( comb( jVictim) ) ];
                suspectInfo.SpotSigmaXYErr = [ suspectInfo.SpotSigmaXYErr, suspectInfo.SpotSigmaXYErr( comb( jVictim) ) ];
                suspectInfo.SpotSigmaZErr = [ suspectInfo.SpotSigmaZErr, suspectInfo.SpotSigmaZErr( comb( jVictim) ) ];
                suspectInfo.num = NParticles;
                
            end
            
            suspectInfoCell{ jComb} = suspectInfo;
            
            % calculate amplitude bias for cost matrix calculation
            biasAmp = zeros( 1, numShadows);
            for jVictim = 1 : numShadows
                
                ampC = suspectInfo.SpotAmp( comb( jVictim) );
                ampO = suspectInfo.SpotAmp; ampO ( comb( jVictim) ) = []; ampO = mean( ampO);
                
                biasAmp( jVictim) = 1 + 0.5 * abs( ampC - ampO ) / ampO; 
                
            end
            
            % Find cost matrix. Use Kalman propagated positions of spots
            % from previous frame
            % find cost matrix from nearest neighbor distance calculation
            
            movieInfo2( refFrame) = refInfo; %6 particles
            movieInfo2( suspectFrame) = suspectInfo; % n < 6 particles
            
            eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
                ' = ' 'costMatRandomDirectedSwitchingMotionLinkSA' '(movieInfo2,kalmanFilterInfo(iFrame),'...
                'costMatParam,probDim,featLifetime,iFrame);'])
            
            nnDist = createDistanceMatrix( suspectInfo.allCoord(:,1:2:end),...
            refInfo.allCoord(:,1:2:end));

            % add bias to nnDist
            for jVictim = 1 : numShadows
            
                nnDist( numConcealers + jVictim, :) = biasAmp( jVictim) .* nnDist( numConcealers + jVictim, :);
                
            end

            % find the effective minimum linking cost for unique matchings
            % Will have to try all combinations of unique matches, for best
            % global match
            pp = perms( 1 : numTrue);
            costs = zeros( size( pp,1), 1);
            for jMatch = 1 : size( pp, 1)
                
                idx2 = pp( jMatch, :);
                idx1 = 1 : length( idx2);
                
                costs( jMatch) = sum( nnDist( sub2ind(size(nnDist), idx1, idx2) ) );
                
            end
            
            % pick the combination with the min overall cost
            [costEff( jComb), idx] = min( costs);
            matches{ jComb} = [ idx1' , [ pp( idx, :)]' ];
                
        end
        
        % We have our first linking assignment between frame1 and frame 2
        [~, idx] = min( costEff);
        match = matches{ idx};
        
        link12 = match(:, 2);
        [~, idxx] = sort( link12, 'ascend');
        link21 = match( idxx,1);
        
        % Store identity of victim in movieInfo2
        movieInfo2( suspectFrame) = suspectInfoCell{ idx};
        movieInfo2( refFrame) = movieInfo( refFrame);
           
        
    end

%     if numFeaturesFrame1 ~= 0 %if there are features in 1st frame
%         
%         if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
%             
%             %calculate cost matrix
%             % -- USER DEFINED FUNCTION -- %
%             eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
%                 ' = ' costMatName '(movieInfo,kalmanFilterInfo(iFrame),'...
%                 'costMatParam,nnDistFeatures(1:numFeaturesFrame1,:),'...
%                 'probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame);'])
% 
%             % % %             %for paper - get number of potential links per feature
%             % % %             numPotLinksPerFeature = [numPotLinksPerFeature; sum(...
%             % % %                 costMat(1:numFeaturesFrame1,1:numFeaturesFrame2)...
%             % % %                 ~=nonlinkMarker,2)];
% 
%             if any(costMat(:)~=nonlinkMarker) %if there are potential links
% 
%                 %link features based on cost matrix, allowing for birth and death
%                 [link12,link21] = lap(costMat,nonlinkMarker,0);
% 
%                 %get indices of features in 2nd frame that are connected to features in 1st frame
%                 indx2C = find(link21(1:numFeaturesFrame2)<=numFeaturesFrame1);
% 
%                 %get indices of corresponding features in 1st frame
%                 indx1C = link21(indx2C);
% 
%                 %find existing tracks that are not connected to features in 2nd frame
%                 numExistTracks = size(trackedFeatureIndx,1);
%                 indx1U = setdiff(1:numExistTracks,indx1C);
%                 numRows = length(indx1U);
%                 
%                 %determine where to store these tracks in auxiliary matrix
%                 %extend auxiliary matrices if necessary
%                 rowStart = rowEnd - numRows + 1;
%                 if rowStart <= 1
%                     trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
%                         trackedFeatureIndxAux];
%                     nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
%                         nnDistFeaturesAux];
%                     prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
%                         prevCostAux];
%                     rowEnd = rowEnd + numTracksWorstCase;
%                     rowStart = rowStart + numTracksWorstCase;
%                 end
%                 
%                 %move rows of tracks that are not connected to points in
%                 %2nd frame to auxilary matrix
%                 trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx(indx1U,:);
%                 
%                 %assign space for new connectivity matrix
%                 tmp = zeros(numFeaturesFrame2,iFrame+1);
% 
%                 %fill in the feature numbers in 2nd frame
%                 tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
% 
%                 %shuffle existing tracks to get the correct connectivity with 2nd frame
%                 tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);
%                 
%                 %update the connectivity matrix "trackedFeatureIndx"
%                 trackedFeatureIndx = tmp;
% 
%                 %repeat for the matrix of nearest neighbor distances
%                 nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures(indx1U,:);
%                 tmp = NaN(numFeaturesFrame2,iFrame+1);
%                 tmp(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;
%                 tmp(indx2C,1:iFrame) = nnDistFeatures(indx1C,:);
%                 nnDistFeatures = tmp;
% 
%                 %repeat for the matrix of linking costs
%                 prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost(indx1U,:);
%                 tmp = NaN(numFeaturesFrame2,iFrame+1);
%                 for i = 1 : length(indx2C)
%                     tmp(indx2C(i),iFrame+1) = costMat(indx1C(i),indx2C(i));
%                 end
%                 tmp(indx2C,1:iFrame) = prevCost(indx1C,:);
%                 prevCost = tmp;
%                 
%                 %update rowEnd to indicate until which row the auxiliary
%                 %matrices are ampty
%                 rowEnd = rowStart - 1;
% 
%                 %get track lifetimes for features in 2nd frame
%                 featLifetime = ones(numFeaturesFrame2,1);
%                 for iFeat = 1 : numFeaturesFrame2
%                     featLifetime(iFeat) = length(find(trackedFeatureIndx(iFeat,:)~=0));
%                 end
% 
%                 %use the Kalman gain from linking to get better estimates
%                 %of the state vector and its covariance matrix in 2nd frame
%                 %as well as state noise and its variance
%                 if selfAdaptive
% 
%                     % -- USER DEFINED FUNCTION -- %
%                     if usePriorInfo %if prior information is supplied
%                         eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
%                             '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
%                             'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
%                             'propagationScheme,kalmanFilterInfo,probDim,'...
%                             'filterInfoPrev(iFrame+1),costMatParam,kalmanFunctions.initialize);'])
%                     else %if no prior information is supplied
%                         eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
%                             '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
%                             'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
%                             'propagationScheme,kalmanFilterInfo,probDim,'...
%                             '[],costMatParam,kalmanFunctions.initialize);'])
%                     end
%                     
%                 end
%                 
%             else %if there are no potential links
%                 
%                 %determine where to store the tracks up to 1st frame in
%                 %auxiliary matrix
%                 %extend auxiliary matrices if necessary
%                 numRows = size(trackedFeatureIndx,1);
%                 rowStart = rowEnd - numRows + 1;
%                 if rowStart <= 1
%                     trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
%                         trackedFeatureIndxAux];
%                     nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
%                         nnDistFeaturesAux];
%                     prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
%                         prevCostAux];
%                     rowEnd = rowEnd + numTracksWorstCase;
%                     rowStart = rowStart + numTracksWorstCase;
%                 end
%                 
%                 %move tracks upto 1st frame to auxiliary matrix
%                 trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;
% 
%                 %assign space for new connectivity matrix
%                 trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
% 
%                 %fill in the feature numbers in 2nd frame
%                 trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
%                 
%                 %repeat for the matrix of nearest neighbor distances
%                 nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures;
%                 nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
%                 nnDistFeatures(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;
% 
%                 %repeat for the matrix of linking costs
%                 prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost;
%                 prevCost = NaN(numFeaturesFrame2,iFrame+1);
%                 
%                 %update rowEnd to indicate until which row the auxiliary
%                 %matrices are ampty
%                 rowEnd = rowStart - 1;
% 
%                 %assign track lifetimes for features in 2nd frame
%                 featLifetime = ones(numFeaturesFrame2,1);
% 
%                 %initialize Kalman filter for features in 2nd frame
%                 if selfAdaptive
% 
%                     % -- USER DEFINED FUNCTION -- %
%                     if usePriorInfo %use a priori information if available
%                         kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
%                         kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
%                         kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
%                     else
%                         eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
%                             '(movieInfo(iFrame+1),probDim,costMatParam);'])
%                         kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
%                         kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
%                         kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
%                     end
%                     
%                 end
% 
%             end %(if any(costMat(:)~=nonlinkMarker))
%             
%         else %if there are no features in 2nd frame
%             
%             %determine where to store the tracks up to 1st frame in
%             %auxiliary matrix
%             %extend auxiliary matrices if necessary
%             numRows = size(trackedFeatureIndx,1);
%             rowStart = rowEnd - numRows + 1;
%             if rowStart <= 1
%                 trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
%                     trackedFeatureIndxAux];
%                 nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
%                     nnDistFeaturesAux];
%                 prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
%                     prevCostAux];
%                 rowEnd = rowEnd + numTracksWorstCase;
%                 rowStart = rowStart + numTracksWorstCase;
%             end
%             
%             %move tracks upto 1st frame to auxiliary matrix
%             trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;
%             
%             %update the connectivity matrix "trackedFeatureIndx"
%             trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
%             
%             %repeat for the matrix of nearest neighbor distances
%             nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures;
%             nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
%             
%             %repeat for the matrix of linking costs
%             prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost;
%             prevCost = NaN(numFeaturesFrame2,iFrame+1);
%             
%             %update rowEnd to indicate until which row the auxiliary
%             %matrices are ampty
%             rowEnd = rowStart - 1;
%             
%             %assign track lifetimes for features in 2nd frame
%             featLifetime = [];
%             
%         end %(if numFeaturesFrame2 ~= 0 ... else ...)
%         
%     else %if there are no features in 1st frame
%         
%         if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
%             
%             %assign space for new connectivity matrix
%             trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
%             
%             %fill in the feature numbers in 2nd frame
%             trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
%             
%             %repeat for the matrix of nearest neighbor distances
%             nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
%             nnDistFeatures(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;
%             
%             %repeat for the matrix of linking costs
%             prevCost = NaN(numFeaturesFrame2,iFrame+1);
%             
%             %assign track lifetimes for features in 2nd frame
%             featLifetime = ones(numFeaturesFrame2,1);
%             
%             %initialize Kalman filter for features in 2nd frame
%             if selfAdaptive
%                 
%                 % -- USER DEFINED FUNCTION -- %
%                 if usePriorInfo %use a priori information if available
%                     kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
%                     kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
%                     kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
%                 else
%                     eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
%                         '(movieInfo(iFrame+1),probDim,costMatParam);'])
%                     kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
%                     kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
%                     kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
%                 end
%                 
%             end
% 
%         else %if there are no features in 2nd frame
% 
%             %assign space for new connectivity matrix
%             trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
%             
%             %repeat for the matrix of nearest neighbor distances
%             nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
%             
%             %repeat for the matrix of linking costs
%             prevCost = NaN(numFeaturesFrame2,iFrame+1);
%             
%             %assign track lifetimes for features in 2nd frame
%             featLifetime = [];
% 
%         end %(if numFeaturesFrame2 ~= 0 ... else ...)
% 
%     end %(if numFeaturesFrame1 ~= 0 ... else ...)
% 
%     %update structure of previous costs
%     prevCostStruct.all = prevCost;
%     prevCostStruct.max = max([prevCostStruct.max; prevCost(:,end)]);
%     prevCostStruct.allAux = prevCostAux;
%     
%     %display progress
%     if verbose
%         progressText(iFrame/(numFrames-1),'Linking frame-to-frame');
%     end
    
end %(for iFrame=1:numFrames-1)

%add information from last frame to auxiliary matrices
numRows = size(trackedFeatureIndx,1);
rowStart = rowEnd - numRows + 1;
if rowStart <= 1
    trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
        trackedFeatureIndxAux];
    nnDistFeaturesAux = [NaN(numRows,numFrames); ...
        nnDistFeaturesAux];
    prevCostAux = [NaN(numRows,numFrames); ...
        prevCostAux];
    rowEnd = rowEnd + numRows;
    rowStart = rowStart + numRows;
end
trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;
nnDistFeaturesAux(rowStart:rowEnd,:) = nnDistFeatures;
prevCostAux(rowStart:rowEnd,:) = prevCost;

%remove all empty rows
trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
clear trackedFeatureIndxAux
nnDistFeatures = nnDistFeaturesAux(rowStart:end,:);
clear nnDistFeaturesAux
prevCost = prevCostAux(rowStart:end,:);
clear prevCostAux

%get total number of tracks
numTracks = size(trackedFeatureIndx,1);

%find the frame where each track begins and then sort the vector
frameStart = zeros(numTracks,1);
for i=1:numTracks
    frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
end
[frameStart,indx] = sort(frameStart);

%rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureIndx = trackedFeatureIndx(indx,:);

%also re-arrange the matrix indicating nearest neighbor distances and
%previouc costs
nnDistFeatures = nnDistFeatures(indx,:);
prevCost = prevCost(indx,:);

%clear some memory
clear costMat tmp

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
%trackedFeatureInfo is in sparse format
trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
    numFrames,probDim);

%take absolute value of all noise variances - this takes care of the
%negative variances used to indicate first appearances
if selfAdaptive && ~usePriorInfo
    for iFrame = 1 : numFrames
        kalmanFilterInfo(iFrame).noiseVar = abs(kalmanFilterInfo(iFrame).noiseVar);
    end
end


%% %%%%% ~~ the end ~~ %%%%%

