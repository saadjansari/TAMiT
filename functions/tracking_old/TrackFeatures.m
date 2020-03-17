function trackResults = TrackFeatures( Movies, movieInfo, trackParams )
% TrackFeatures :  Uses kalman filters to generate feature tracks given
% detected feautures for every frame

    % Specify for which channels to track the detected features.
    trackChannels = [2];
    
    % Models used for channels
    ChannelDetection(1).Model = [3];
    ChannelDetection(2).Model = [1];
    
    % For each selected channel, track the detected features
    for jChannel = trackChannels
         
        % Load the detection data
        fitResults = movieInfo.fitResults( ChannelDetection( jChannel).Model ).frame;
        
        for jFrame = 1 : length(fitResults)
            fitResults(jFrame).amp = [ [ fitResults(jFrame).SpotAmp]' , [fitResults(jFrame).SpotAmpErr]'];
            fitResults(jFrame).zCoord = fitResults(jFrame).SpotZ;
        end

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
    
end

