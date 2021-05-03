classdef TrackLines < TrackFeatures
    % Class for tracking line features in cells.
    % MONOPOLAR CELL
    % Uses the Danuser u-track software and the generic TrackFeatures class
    properties
        startPt % nDim x nTime
        scaling = [0.5,1,3] % [phi,theta,len]
    end
    
    methods (Access = public)
        
        function obj = TrackLines(data_xyzt, times, timestep, sizeVoxels)
            
            obj = obj@TrackFeatures('Line', data_xyzt, times, timestep, sizeVoxels);
            
        end
        
        function obj = parseMainFeature( obj, mainFeatures)
            
            obj.startPt = zeros(3, length( obj.times),2 );
            for jTime = 1 : length( obj.times)
                feat = mainFeatures{jTime};

                % Create a structure with fields xCoord,yCoord,zCoord,amp
%                 % The coordinates will be of the endposition
%                 xC = []; yC = []; zC = []; amp = [];
%                 faster = feat.featureList{1};
%                 obj.startPt(:, jTime, 1) = faster.featureList{1}.position;
%                 obj.startPt(:, jTime, 2) = faster.featureList{1}.err_position;
%                 for jc = 2: faster.numFeatures
%                     cc = faster.featureList{jc}.endPosition;
%                     ccErr = faster.featureList{jc}.err_endPosition;
%                     xC = [xC; [cc(1), ccErr(1)]];
%                     yC = [yC; [cc(2), ccErr(2)]];
%                     zC = [zC; [cc(3), ccErr(3)]];
%                     amp = [amp; [faster.featureList{jc}.amplitude, faster.featureList{jc}.err_amplitude]];
%                 end
%                 movieInfo(jTime).xCoord = xC;
%                 movieInfo(jTime).yCoord = yC;
%                 movieInfo(jTime).zCoord = zC;
%                 movieInfo(jTime).amp = amp;

                % The coordinates will be (L,phi,theta) of the end
                % position.
                xC = []; yC = []; zC = []; amp = [];
                faster = feat.featureList{1};
                obj.startPt(:, jTime, 1) = faster.featureList{1}.position;
                obj.startPt(:, jTime, 2) = faster.featureList{1}.err_position;
                for jc = 2: faster.numFeatures
                    
                    ctheta1 = faster.featureList{jc}.theta(1)/ obj.scaling(1);
                    ctheta2 = faster.featureList{jc}.theta(2)/ obj.scaling(2);
                    clen = faster.featureList{jc}.length / obj.scaling(3);
                    
                    xC = [xC; [ctheta1, 0]];
                    yC = [yC; [ctheta2, 0]];
                    zC = [zC; [clen, 0]];
                    
                    amp = [amp; [faster.featureList{jc}.amplitude, faster.featureList{jc}.err_amplitude]];
                end
                movieInfo(jTime).xCoord = xC;
                movieInfo(jTime).yCoord = yC;
                movieInfo(jTime).zCoord = zC;
                movieInfo(jTime).amp = amp;
            end
            obj.movieInfo = movieInfo;
            
        end
        
        function obj = trackUTRACK( obj)
           
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

            % additional input

            %saveResults
            saveResults.dir = '~/Documents/Projects/ImageAnalysis/SingleCell'; %directory where to save input and output
            saveResults.filename = 'tracksTest.mat'; %name of file where input and output are saved
            % saveResults = 0; %don't save results

            %verbose state
            verbose = 1;

            %problem dimension
            probDim = 3;

            % tracking function call

            [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(obj.movieInfo,...
                costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
            
            obj.tracksFinal = tracksFinal;
        end
        
        function dyfeats = mat2dyfeats( obj, xC, yC, zC, aC)
            % Get dynamic features from matrix xt, yt,zt data
            
            numT = size(xC,2);
            numFeat = size(xC,1);

            % Determine number of distinguisable colors needed
            nFeatsPerFrame = sum( ~isnan(xC(:,:,1)), 1);
            num_col = max( nFeatsPerFrame );
            num_col = numFeat;
            cols = distinguishable_colors(num_col, {'w','k'});
            
            % Get feature start and end times
            not_nans = ~isnan( xC(:,:,1));
            startend = zeros(numFeat,2);
            for jf = 1: numFeat
                startend(jf,1) = find(not_nans(jf,:),1,'first');
                startend(jf,2) = find(not_nans(jf,:),1,'last');
            end
            feat_timespan = startend(:,2)-startend(:,1)+1;
            
            % Initialize dynamic features
            dyfeats = cell( 1, numFeat);
            for jf = 1: numFeat
                
                % start and end times
                t0 = startend(jf,1);
                t1 = startend(jf,2);
                
                % get amplitudes
                amp = aC( jf, t0:t1, :);
                
%                 % get position matrix
%                 pos = zeros(2,3, feat_timespan(jf) ,2 );
%                 pos(1,:,:,:) = obj.startPt(:,t0:t1,:);
%                 pos(2,1,:,:) = xC( jf,t0:t1,:);
%                 pos(2,2,:,:) = yC( jf,t0:t1,:);
%                 pos(2,3,:,:) = zC( jf,t0:t1,:);
%                 
%                 dyfeats{jf} = DynamicLine( startend(jf,1), startend(jf,2), pos, amp, cols( 1+rem(jf,num_col), :), obj.timeStep, obj.sizeVoxels ); 

                
                % get position matrix
                pos0 = obj.startPt(:,t0:t1,:);
                
                theta = zeros(2, feat_timespan(jf), 2);
                theta(1,:,:) = xC( jf,t0:t1,:)*obj.scaling(1);
                theta(2,:,:) = yC( jf,t0:t1,:)*obj.scaling(2);
                lens = zC( jf,t0:t1,:)*obj.scaling(3);
                
                dyfeats{jf} = DynamicLine_sph( startend(jf,1), startend(jf,2), pos0, lens, theta, amp, cols( 1+rem(jf,num_col), :), obj.timeStep, obj.sizeVoxels ); 
            
            end
            
        end
        
    end
    
    
end