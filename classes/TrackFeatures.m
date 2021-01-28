classdef TrackFeatures
    % Class for tracking features in cells
    % Uses the Danuser u-track software
    properties
        type
        data_xyt
        movieInfo
        tracksFinal
        times = []
        features
        timeStep
        sizeVoxels = [0.1067 0.1067 0.5]; %default
        numVoxels 
        playMovie = 1
        saveMovie = 0
    end
    
    methods (Access = public)
        
        function obj = TrackFeatures(type, data_xyt, times, timeStep, sizeVoxels)
            
            if ~strcmp(type, 'Line') && ~strcmp( type, 'Curve') && ~strcmp( type, 'Spot')
                error('TrackFeature: unknown type, allowed values are Line, Curve, Spot')
            end

            obj.type = type;
            obj.data_xyt = data_xyt;
            obj.times = times;
            obj.timeStep = timeStep;
            obj.sizeVoxels = sizeVoxels; 
            obj.numVoxels = size(obj.data_xyt);
            addpath( genpath('functions/external/u-track') )
            
        end
        
        function obj = parseTracksFinal( obj)
            % Parse tracking results and create dynamic features
            
            
            numFrames = size(obj.data_xyt,3);  
            minLength = 3;
            
            % store track positions, get track status and point status

            %get number of tracks
            numTracks = length(obj.tracksFinal);
            %get track start and end times
            trackSEL = getTrackSEL(obj.tracksFinal);

            %get number of segments making each track
            numSegments = zeros(numTracks,1);
            for i = 1 : numTracks
                numSegments(i) = size(obj.tracksFinal(i).tracksCoordAmpCG,1);
            end

            %locate the row of the first segment of each compound track in the
            %big matrices of all tracks (to be constructed in the next step)
            trackStartRow = ones(numTracks,1);
            for iTrack = 2 : numTracks
                trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
            end

            %find total number of segments in all tracks (i.e. number of rows in big
            %matrices)
            numSegmentsTracks = trackStartRow(end)+numSegments(end)-1;

            %put all tracks together in one big matrix
            %put the x-coordinates in one matrix and the y-coordinates in another
            %indicate the status of each point
            xCoordMatAll = NaN*ones(numSegmentsTracks,numFrames,2);
            yCoordMatAll = xCoordMatAll;
            zCoordMatAll = xCoordMatAll;
            ampMatAll = xCoordMatAll;
            
            xCoordMatAll_err = NaN*ones(numSegmentsTracks,numFrames);
            yCoordMatAll_err = xCoordMatAll_err;
            zCoordMatAll_err = xCoordMatAll_err;
            
            for iTrack = 1 : numTracks

                %get track start and end times
                startTime = trackSEL(iTrack,1);
                endTime   = trackSEL(iTrack,2);

                %store x-coordinates + errors
                xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,1) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
                xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,2) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,5:8:end);

                %store y-coordinates
                yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,1) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
                yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,2) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,6:8:end);
                
                %store z-coordinates
                zCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,1) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,3:8:end);
                zCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,2) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,7:8:end);
                
                %store amplitudes
                ampMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,1) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,4:8:end);
                ampMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
                    numSegments(iTrack)-1,startTime:endTime,2) = ...
                    obj.tracksFinal(iTrack).tracksCoordAmpCG(:,8:8:end);
                
                

                %get sequence of events of track
                seqOfEvents = obj.tracksFinal(iTrack).seqOfEvents;

                %assign point status for features just after and before a split
                %also, in the frame just before splitting, give the
                %splitting track the position of the track it split from
                points2consider = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)))';
                for iPoint = points2consider
                    xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1)-1 ,:) = xCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1 , :);
                    yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1)-1 , :) = yCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1 , :);
                    zCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1)-1, :) = zCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1, :);

                end

                %assign point status for features just before and after a merge
                %also, in the frame just after merging, give the
                %merging track the position of the track it merged from
                points2consider = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)))';
                for iPoint = points2consider
                    
                    xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1), :) = xCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1), :);
                    yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1), :) = yCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1), :);
                    zCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                        seqOfEvents(iPoint,1), :) = zCoordMatAll(trackStartRow(iTrack)...
                        +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1), :);
                end

            end %(for iTrack = 1 : numTracks)

            % Call specialized function to get dynamic features
            obj.features = obj.mat2dyfeats( xCoordMatAll, yCoordMatAll, zCoordMatAll, ampMatAll);
            
            idxRm = [];
            for idx = 1 : length(obj.features)
                if obj.features{idx}.time_end-obj.features{idx}.time_start <= minLength
                    idxRm =[idxRm, idx];
                end
            end
            obj.features(idxRm) = [];
            
        end
        
        
    end
    
    methods (Abstract = true)
        
        obj = mat2dyfeats( obj, xC, yC, zC)

    end
    
    
end