function F = plotTracks( imData4, fitResults, tracksFinal)
% plotTracks : overlay the features from trackResults on MIP images from
% imData4
% Assigns a unique color to each feature in trackResults

    % Lets start by finding the start and end frames for which we have
    % features
    startFrame = [];
    endFrame = [];
    for jFrame = 1 : length(fitResults)
        
        if ~isempty( fitResults(jFrame).Bkg) && isempty( startFrame)
            startFrame = jFrame;
        end
        if ~isempty( fitResults(jFrame).Bkg)
            endFrame = jFrame;
        end
        
    end
%     startFrame = find( [fitResults.Bkg] >= 0, 1, 'first');
%     endFrame = find( [fitResults.Bkg] >= 0 , 1, 'last');
    numFrames = endFrame - startFrame + 1;

    % Now lets make the MIP images of each of these frames and stroe them
    % in a cell array
    imMIPCell = cell( 1, numFrames);
    for jFrame = 1 : numFrames
        
        currFrame = startFrame + jFrame - 1;
        imMIPCell{ jFrame} = max( imData4( :, :, :, currFrame), [], 3);
        
    end
    
    % Now lets process the features
    % Extract number of unique features (which may or may not last the
    % entire duration of the movie)
    nFeatures = length( tracksFinal );
    % Give a unique color to each feature
    cFeatures = distinguishable_colors( nFeatures, [ 0, 0, 0]);
    
    % For each feature, determine what frames it exists in, so we'll only
    % plot it in those frames
    for jFeature = 1 : nFeatures
        
       % find number of tracks in this compound track
       ntracks = size( tracksFinal( jFeature).tracksFeatIndxCG, 1);
       
       % find start and end of each track
       framesTrack = [];
       for jT = 1 : ntracks
           % find rows corresponding to this track in seqOfEvents
           idx = find( tracksFinal( jFeature).seqOfEvents(:, 3) == jT, 2);
           startTrack = tracksFinal( jFeature).seqOfEvents( idx(1), 1);
           endTrack = tracksFinal( jFeature).seqOfEvents( idx(2), 1);
           framesTrack = [ framesTrack, startTrack : endTrack ];
       end
       
       frames = unique( framesTrack);
       frames = frames(1) : frames(end);
       FeatureData( jFeature).frames = frames;
       
       % Now we'll record the (x,y) locations of this feature for each of
       % the defined frames
       for jFrame = 1 : length( frames)
           
           currFrame = startFrame + frames(1) + jFrame - 2;
           featNum = tracksFinal( jFeature).tracksFeatIndxCG( :, jFrame);
           featNum( featNum == 0) = [];
           featNum( isnan(featNum)) = [];
           if ~isempty( featNum)
               FeatureData( jFeature).x{ jFrame} = fitResults( currFrame).SpotX( featNum);
               FeatureData( jFeature).y{ jFrame} = fitResults( currFrame).SpotY( featNum);
           end
       end
       
    end
    
    % Now lets parse this data.
    for jFrame = 1 : numFrames
        
        FrameData( jFrame).image = imMIPCell{ jFrame};
        FrameData( jFrame).x = [];
        FrameData( jFrame).y = [];
        FrameData( jFrame).colors = [];
        % Check if any features exist in this frame and record their (x,y)
        % positions
        feats = [];
        for jFeat = 1 : nFeatures
            if any( FeatureData( jFeat).frames == jFrame)
                
                idx = find( FeatureData( jFeat).frames == jFrame );
                nfeats = length( FeatureData( jFeat).x{idx});
                FrameData( jFrame).x = [ FrameData( jFrame).x, FeatureData( jFeat).x{idx} ];
                FrameData( jFrame).y = [ FrameData( jFrame).y, FeatureData( jFeat).y{idx} ];
                FrameData( jFrame).colors = [ FrameData( jFrame).colors; repmat( cFeatures( jFeat, :), nfeats, 1) ];

            end
        end
    end
    
    % Now we can go ahead with the plotting and movie-making
    
    F(numFrames) = struct('cdata',[],'colormap',[]);
    
    for jFrame = 1 : numFrames
        h = figure;
        hold on
        imagesc( FrameData( jFrame).image); colormap gray; axis equal
        
        for jPlot = 1 : length( FrameData( jFrame).x )
            plot( FrameData( jFrame).x( jPlot), FrameData( jFrame).y( jPlot), ...
                'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, ...
                'color', FrameData( jFrame).colors( jPlot, :) );
        end
        hold off
        title( ['Tracked Kinetochores : Frame ', num2str(jFrame) ])
        F( jFrame) = getframe( h);
        close(h)
        
    end

end