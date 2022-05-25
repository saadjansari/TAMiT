classdef TrackSpots < TrackFeatures
    % Class for tracking spot features in cells.
    % MITOTIC CELL
    % Uses the Danuser u-track software and the generic TrackFeatures class
    properties
        scaling = [3,3,1] % [position]
    end
    
    methods (Access = public)
        
        function obj = TrackSpots(data_xyzt, times, timestep, sizeVoxels)
            
            obj = obj@TrackFeatures('Spot', data_xyzt, times, timestep, sizeVoxels);
            
        end
        
        function obj = parseMainFeature( obj, mainFeatures)
            
            for jTime = 1 : length( obj.times)
                feat = mainFeatures{jTime};

                % The coordinates will be position.
                xC = []; yC = []; zC = []; amp = [];
                for jc = 1: feat.numFeatures
                    
                    xx = feat.featureList{jc}.position(1)/ obj.scaling(1);
                    yy = feat.featureList{jc}.position(2)/ obj.scaling(2);
                    zz = feat.featureList{jc}.position(3) / obj.scaling(3);
                    
                    xC = [xC; [xx, 0]];
                    yC = [yC; [yy, 0]];
                    zC = [zC; [zz, 0]];
                    
                    amp = [amp; [feat.featureList{jc}.amplitude, feat.featureList{jc}.err_amplitude]];
                end
                movieInfo(jTime).xCoord = xC;
                movieInfo(jTime).yCoord = yC;
                movieInfo(jTime).zCoord = zC;
                movieInfo(jTime).amp = amp;
            end
            obj.movieInfo = movieInfo;
            
        end
        
        function obj = trackUTRACK( obj)
           
            gapCloseParam.timeWindow = 10; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

            %optional input:
            gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            %function name
            costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

            %parameters

            parameters.linearMotion = 2; %use linear motion Kalman filter.
            parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            parameters.maxSearchRadius = 5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            parameters.brownStdMult = 5; %multiplication factor to calculate search radius from standard deviation.

            parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
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
%             saveResults.dir = '~/Documents/Projects/ImageAnalysis/SingleCell'; %directory where to save input and output
%             saveResults.filename = 'tracksTest.mat'; %name of file where input and output are saved
            saveResults = 0; %don't save results

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
                pos0 = zeros(3, feat_timespan(jf), 2);
                
                pos0(1,:,:) = xC( jf,t0:t1,:)*obj.scaling(1);
                pos0(2,:,:) = yC( jf,t0:t1,:)*obj.scaling(2);
                pos0(3,:,:) = zC( jf,t0:t1,:)*obj.scaling(3);
                
                dyfeats{jf} = DynamicSpot( startend(jf,1), startend(jf,2), pos0, amp, cols( 1+rem(jf,num_col), :), obj.timeStep, obj.sizeVoxels ); 
            
            end
            
        end
        
    end
    
    
end