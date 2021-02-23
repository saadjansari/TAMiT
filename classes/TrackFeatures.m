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
        
        function obj = parseTracksFinal( obj, feats)
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
            
            % Assign feature objects from detection process, to tracked
            % features. This enables filling in parameters that werent used
            % in the tracking process, and helps with things like
            % displaying figures and making movies
            obj = obj.matchFeatures( feats);
            
        end
        
        function obj = matchFeatures( obj, feats)
            % Matches features from detection to features from tracking
            
            for jf = 1 : length( obj.features)
                cf = obj.features{jf};
                obj.features{jf} = cf.matchFeatures( feats(cf.time_start : cf.time_end) );
            end
                
        end
        
        function drawMatchedFeature( obj, ax, jtime)
            % Draw matched detected feature
            
            for jf = 1 : length(obj.features)
            
                if obj.features{jf}.existAtTime(jtime)
                    obj.features{jf}.matched_feats{ 1+ jtime - obj.features{jf}.time_start}.displayFeature( ax);
                end
            end
            
        end
        
        function [lbls, cols] = drawMatchedFeature_withlabel( obj, ax, jtime, lbls, cols, offset)
            % Draw matched detected feature along with its label

            for jf = 1 : length(obj.features)
            
                if obj.features{jf}.existAtTime(jtime)
                    cf = obj.features{jf}.matched_feats{ 1+ jtime - obj.features{jf}.time_start};
                    if ~isempty( cf)
                        cf.displayFeature( ax);
                        % feature label and color
                        lbls = { lbls{:}, num2str(jf+offset) };
                        cols = [ cols ; cf.display{2} ];
                    end
                    
                end
            end
            
        end
        
    end
    
    methods( Access = public, Static = true)
       
        function makeMovie( type, img, times, feats, movpath)
            % Make a tracking movie
            
            do_save_movie = 1;
            draw_feat_label = 1;
            
            if do_save_movie
                mov = cell(1,length(times) );
            end
            
            
            % Make a figure
            f = figure(); f.Color = 'white';
            h = tight_subplot(1,2, 0.001, 0.125,0.11);
            
            for jt = 1 : length(times)
                
                if jt ~= 1
                    cla( h(1) ); % clear axes if not first frame
                    cla( h(2) ); % clear axes if not first frame
                    legend off
                end
                
                jax = 1;
                set(f, 'currentaxes', h(jax) );
                % Display image
                imagesc( img(:,:,jt) ); colormap gray; axis equal; 
                xlim( [1 size(img,2) ]); ylim( [1 size(img,1) ]); 
                set( h(jax), 'Xtick', [], 'Ytick', []); hold on;
                
                jax = 2;
                set(f, 'currentaxes', h(jax) );
                % Display image
                imagesc( img(:,:,jt) ); colormap gray; axis equal; 
                xlim( [1 size(img,2) ]); ylim( [1 size(img,1) ]); 
                set( h(jax), 'Xtick', [], 'Ytick', []); hold on;

                % Display features
                switch type
                    case 'MitosisBud'

                        % Display the spindles
                        spindles = feats{1};
                        spindles{jt}.displayFeature( h(jax));
                        
                        bud1 = feats{2};
                        bud2 = feats{3};
                        
                        if draw_feat_label
                            % Display features at each SPB
                            [lbl, cols] = bud1.drawMatchedFeature_withlabel( h(jax), jt, {}, [], 0);                            
                            [lbl, cols] = bud2.drawMatchedFeature_withlabel( h(jax), jt, lbl, cols, length( bud1.features) );
                            % draw feat labels
                            if ~isempty(lbl)
                                clear p
                                for ii = 1:size(cols,1)
                                    p(ii) = patch(NaN, NaN, cols(ii,:));
                                end
                                try
                                legend(p, lbl);
                                catch
                                    stopph=1;
                                end
                            end
                        else
                            % Display features at each SPB
                            bud1.drawMatchedFeature( h(jax), jt);                            
                            bud2.drawMatchedFeature( h(jax), jt);
                        end
                        

                    case 'Mitosis'
                        error('not set up yet')
                        % Display the spindles
                        spindles = feats{1};
                        spindles{jt}.displayFeature( h(jax))

                    case 'Monopolar'
                        
                        % Display the SPBs
                        spbs = feats{1};
                        spbs{jt}.displayFeature( h(jax));
                        
                        % Display features at 1st SPB
                        lines = feats{2};
                        lines.drawMatchedFeature( h(jax), jt);
                        

                end
                title(['Time = ', num2str( times(jt))]);
                drawnow;
                pause(0.1);
                
                if do_save_movie
                    mov{ jt} = getframe(f);
                end
            end
            
            if do_save_movie
                mname = 'tracks';
                writerObj = VideoWriter([movpath, filesep, mname],'MPEG-4');
                writerObj.FrameRate = 5;
                writerObj.Quality = 100;

                % open the video writer
                open(writerObj);

                % write the frames to the video
                for frame = 1 : length( mov)
                    % convert the image to a frame
                    writeVideo( writerObj, mov{frame} );
                end

                % close the writer object
                close(writerObj);
            end
            close(f);
            
        end
        
    end
    
    methods (Abstract = true)
        
        % Converts data from matrices to dynamic feature objects
        obj = mat2dyfeats( obj, xC, yC, zC)
        
    end
    
    
end