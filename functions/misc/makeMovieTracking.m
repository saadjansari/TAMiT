function makeMovieTracking(tracksFinal, movieMat, startPt, ...
    saveMovie,movieName,colorTracks,dir2saveMovie,minLength,...
    plotFullScreen, movieType)
%OVERLAYTRACKSMOVIENEW overlays tracks obtained via trackCloseGapsKalman on movies with variable color-coding schemes
%
%SYNPOSIS overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
%    saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
%    imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
%    colorTracks,firstImageFile,dir2saveMovieminLength,plotFullScreen,...
%    movieType)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe].
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]
%       dragtailLength: Length of drag tail (in frames).
%                       Optional. Default: 10 frames.
%                       ** If dragtailLength = 0, then no dragtail.
%                       ** To show tracks from their beginning to their end,
%                       set dragtailLength to any value longer than the
%                       movie.
%                       ** To show tracks statically while features dance
%                       on them, use -1.
%                       ** To show tracks from their beginning to their
%                       end, and to retain tracks even after the particle
%                       disappears, use -2.
%       saveMovie     : 1 to save movie (as Quicktime), 0 otherwise.
%                       Optional. Default: 0.
%       movieName     : filename for saving movie.
%                       Optional. Default: TrackMovie (if saveMovie = 1).
%       filterSigma   : 0 to overlay on raw image, PSF sigma to overlay on
%                       image filtered with given filterSigma.
%                       Optional. Default: 0.
%       classifyGaps  : 1 to classify gaps as "good" and "bad", depending
%                       on their length relative to the legnths of the
%                       segments they connect, 0 otherwise.
%                       Optional. Default: 1.
%       highlightES   : 1 to highlight track ends and starts, 0 otherwise.
%                       Optional. Default: 1.
%       showRaw       : 1 to add raw movie to the left of the movie with
%                       tracks overlaid, 2 to add raw movie at the top of
%                       the movie with tracks overlaid, 0 otherwise.
%                       Optional. Default: 0.
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%       onlyTracks    : 1 to show only tracks without any symbols showing
%                       detections, closed gaps, merges and splits; 0 to
%                       show symbols on top of tracks.
%                       Optional. Default: 0.
%       classifyLft   : 1 to classify objects based on (1) whether they
%                       exist throughout the whole movie, (2) whether they
%                       appear OR disappear, and (3) whether they appear
%                       AND disappear; 0 otherwise.
%                       Optional. Default: 1.
%       diffAnalysisRes:Diffusion analysis results (either output of
%                       trackDiffusionAnalysis1 or trackTransientDiffusionAnalysis2).
%                       Needed if tracks/track segments are to be
%                       colored based on their diffusion classification.
%                       With this option, classifyGaps, highlightES and
%                       classifyLft are force-set to zero, regardless of input.
%                       Optional. Default: None.
%       intensityScale: 0 to autoscale every image in the movie, 1
%                       to have a fixed scale using intensity mean and std,
%                       2 to have a fixed scale using minimum and maximum
%                       intensities.
%                       Optional. Default: 1.
%       colorTracks   : 1 to color tracks by rotating through 7 different
%                       colors, 0 otherwise. With this option,
%                       classifyGaps, highlightES and classifyLft are
%                       force-set to zero, regardless of input.
%                       Option ignored if diffAnalysisRes is supplied.
%                       Optional. Default: 0.
%       movieMat      : Specify matrix containing movie (XYT)
%       dir2saveMovie:  Directory where to save output movie.
%                       If not input, movie will be saved in directory where
%                       images are located.
%                       Optional. Default: [].
%       minLength     : Minimum length of tracks to be ploted.
%                       Optional. Default: 1.
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%       movieType     : 'mov' to make a Quicktime movie using MakeQTMovie,
%                       'avi' to make AVI movie using Matlab's movie2avi,
%                       'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                       using ImageMagick and ffmpeg. These options work
%                       only under linux or mac.
%                       Optional. Default: 'mov'.
%
%OUTPUT the movie.
%
%REMARKS Color-coding:
%        ** Without diffusion classification, all tracks have a neutral
%        color, while objects are color coded in the following way:
%               * Detected object just after appearance: Green circle.
%               * Detected object just before disappearance: Yellow
%                 circle.
%               * Detected object in middle of trajectory that spans
%                 whole movie: White circle.
%               * Detected object in middle of trajectory that appears OR
%                 disappears within movie: Magenta circle.
%               * Detected object in middle of trajectory that appears AND
%                 disappears within movie: Red circle.
%               * Gap that is short than both segments it connects: Cyan
%                 star.
%               * Gap that is longer than at least one ofthe segments it
%                 connects: Blue star.
%               * Object before and after splitting: Green diamond.
%               * OBject before and after merging: Yellow diamond.
%           When classifyGaps = 0, all gaps are cyan.
%           When highlighES = 0, no green and yellow circles.
%           When classifyLft = 0, all objets in middle of trajectory are white.
%
%       ** With diffusion classification, all objects and gaps have neutral
%       color (merges and splits are diamonds), while tracks and track
%       segments are color-coded in the following way:
%               * Type 1: Linear + 1D confined diffusion: Orange.
%               * Type 2: Linear + 1D normal diffusion: Red.
%               * Type 3: Linear + 1D super diffusion: Green.
%               * Type 4: Linear + too short for 1D classification: Yellow.
%               * Type 5: Random/Unclassified + 2D confined diffusion: Blue.
%               * Type 6: Random/Unclassified + 2D normal diffusion: Cyan.
%               * Type 7: Random/Unclassified + 2D super diffusion: Magenta.
%               * Type 8: Random + too short for 2D classification: Purple.
%               * Type 0: Too short for any analysis: Light pink.
%
%Khuloud Jaqaman, August 2007
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% input - basic

%check whether correct number of input arguments was used
if nargin < 1
    disp('--overlayTracksMovieNew: Incorrect number of input arguments!');
    return
end

numFrames = size(movieMat,3);  
isx = size(movieMat,1);
isy = size(movieMat,2);

%keep only tracks with minimum requested length
if isempty(minLength)
    minLength = 1;
end
if minLength > 1
    criteria.lifeTime.min = minLength;
    indx = chooseTracks(tracksFinal,criteria);
    tracksFinal = tracksFinal(indx,:);
end

%get first and last frames where there are tracks
allEvents = vertcat(tracksFinal.seqOfEvents);
tracksFirstFrame = min(allEvents(:,1));
tracksLastFrame = max(allEvents(:,1));

%check startend and assign default if necessary
startend = [tracksFirstFrame tracksLastFrame];

%get number of frames in movie to be made
numFramesMovie = diff(startend) + 1;

%check whether an area of interest was input
imageRange = [1 isx; 1 isy];

%% input - additional parameters

%check whether to save movie
if isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && isempty(movieName)
    movieName = 'trackMovie.mov';
end

%check where to save resulting movie
if saveMovie && isempty(dir2saveMovie)
    dir2saveMovie = pwd;
end

%check whether to use full screen for plotting
if isempty(plotFullScreen)
    plotFullScreen = 0;
end

%decide on movie type
if isempty(movieType)
    movieType = 'mov';
end

%define colors to loop through
colorLoop = [1 0.7 0.7; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: 'light pink',r,g,b,y,m,c

%% store track positions, get track status and point status

%get number of tracks
numTracks = length(tracksFinal);

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%give tracks status based on the frames they span:
%2: track exists throughout movie
%1: track exists either in first frame or in last frame
%0: track does not exist in both first frame and last frame
trackStatus  = (trackSEL(:,1) == tracksFirstFrame) + (trackSEL(:,2) == tracksLastFrame);
trackStatus(:) = 0;

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
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

%construct a matrix indicating point status in big matrices:
%-2: bad gap (gap length > either segment length on its sides)
%-1: good gap (gap length < both segment lengths on its sides)
%0 : before track start or after track end
%1 : detected feature in the middle of a track of trackStatus = 0
%2 : detected feature in the middle of a track of trackStatus = 1
%3 : detected feature in the middle of a track of trackStatus = 2
%4 : detected feature just after a birth
%5 : detected feature just before a death
%6 : detected feature just after a split
%7 : detected feature just before a merge
pointStatus = zeros(numSegmentsTracks,numFrames);

%put all tracks together in one big matrix
%put the x-coordinates in one matrix and the y-coordinates in another
%indicate the status of each point
xCoordMatAll = NaN*ones(numSegmentsTracks,numFrames);
yCoordMatAll = xCoordMatAll;
for iTrack = 1 : numTracks
    
    %get track start and end times
    startTime = trackSEL(iTrack,1);
    endTime   = trackSEL(iTrack,2);
    
    %store x-coordinates
    xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    
    %store y-coordinates
    yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    
    %assign point status for features in the middle of the track
    pointStatus(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = trackStatus(iTrack) + 1;
    
    %get sequence of events of track
    seqOfEvents = tracksFinal(iTrack).seqOfEvents;
    
    %assign point status for features just after and before a split
    %also, in the frame just before splitting, give the
    %splitting track the position of the track it split from
    points2consider = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
    end
    
    %assign point status for features just before and after a merge
    %also, in the frame just after merging, give the
    %merging track the position of the track it merged from
    points2consider = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
    end
    
end %(for iTrack = 1 : numTracks)


%% divide tracks based on diffusion analysis or just into groups to be colored separately

%if tracks are to be individually colored
if colorTracks
    
    %divide tracks among 9 matrices that will get their own colors
    
    %x-coordinates ...
    [xCoordMatAll0,xCoordMatAll1,xCoordMatAll2,xCoordMatAll3,xCoordMatAll4,...
        xCoordMatAll5,xCoordMatAll6,xCoordMatAll7,xCoordMatAll8] = deal(NaN(size(xCoordMatAll)));
    xCoordMatAll0(1:9:end,:) = xCoordMatAll(1:9:end,:);
    xCoordMatAll1(2:9:end,:) = xCoordMatAll(2:9:end,:);
    xCoordMatAll2(3:9:end,:) = xCoordMatAll(3:9:end,:);
    xCoordMatAll3(4:9:end,:) = xCoordMatAll(4:9:end,:);
    xCoordMatAll4(5:9:end,:) = xCoordMatAll(5:9:end,:);
    xCoordMatAll5(6:9:end,:) = xCoordMatAll(6:9:end,:);
    xCoordMatAll6(7:9:end,:) = xCoordMatAll(7:9:end,:);
    xCoordMatAll7(8:9:end,:) = xCoordMatAll(8:9:end,:);
    xCoordMatAll8(9:9:end,:) = xCoordMatAll(9:9:end,:);
    
    %y-coordinates ...
    [yCoordMatAll0,yCoordMatAll1,yCoordMatAll2,yCoordMatAll3,yCoordMatAll4,...
        yCoordMatAll5,yCoordMatAll6,yCoordMatAll7,yCoordMatAll8] = deal(NaN(size(yCoordMatAll)));
    yCoordMatAll0(1:9:end,:) = yCoordMatAll(1:9:end,:);
    yCoordMatAll1(2:9:end,:) = yCoordMatAll(2:9:end,:);
    yCoordMatAll2(3:9:end,:) = yCoordMatAll(3:9:end,:);
    yCoordMatAll3(4:9:end,:) = yCoordMatAll(4:9:end,:);
    yCoordMatAll4(5:9:end,:) = yCoordMatAll(5:9:end,:);
    yCoordMatAll5(6:9:end,:) = yCoordMatAll(6:9:end,:);
    yCoordMatAll6(7:9:end,:) = yCoordMatAll(7:9:end,:);
    yCoordMatAll7(8:9:end,:) = yCoordMatAll(8:9:end,:);
    yCoordMatAll8(9:9:end,:) = yCoordMatAll(9:9:end,:);
    
end

%% make movie

%initialize movie if it is to be saved
if saveMovie
    movieVar = struct('cdata',[],'colormap',[]);
    movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%go over all specified frames and find minimum and maximum intensity in all
%of them combined
intensityMinMax = [];

%go over all specified frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h     = figure();
    set(h,'Position',scrsz);
else
    figure
end
for iFrame = 1 : numFramesMovie
            
    %read specified image
    imageStack = movieMat(:,:,iFrame);
    
    %plot image in current frame and show frame number
    clf;
    
    axes('Position',[0 0 0.495 1]);
    imshow(imageStack,intensityMinMax);
    hold on;

    textDeltaCoord = min(diff(imageRange,[],2))/20;

    text(textDeltaCoord,...
        textDeltaCoord,[num2str(((iFrame+startend(1)-1)-1)*0.1,'%4.1f') ' s'],'Color','yellow');
    axes('Position',[0.505 0 0.495 1]);
    imshow(imageStack,intensityMinMax);

    hold on;
            
    %get tracks to plot
    plotOrNot = 0;
    if iFrame > 1
        dragTailStart = 1;
        indx2keep = find(pointStatus(:,iFrame)~=0);
        xCoord2plot0 = (xCoordMatAll0(indx2keep,dragTailStart:iFrame))';
        yCoord2plot0 = (yCoordMatAll0(indx2keep,dragTailStart:iFrame))';
        xCoord2plot1 = (xCoordMatAll1(indx2keep,dragTailStart:iFrame))';
        yCoord2plot1 = (yCoordMatAll1(indx2keep,dragTailStart:iFrame))';
        xCoord2plot2 = (xCoordMatAll2(indx2keep,dragTailStart:iFrame))';
        yCoord2plot2 = (yCoordMatAll2(indx2keep,dragTailStart:iFrame))';
        xCoord2plot3 = (xCoordMatAll3(indx2keep,dragTailStart:iFrame))';
        yCoord2plot3 = (yCoordMatAll3(indx2keep,dragTailStart:iFrame))';
        xCoord2plot4 = (xCoordMatAll4(indx2keep,dragTailStart:iFrame))';
        yCoord2plot4 = (yCoordMatAll4(indx2keep,dragTailStart:iFrame))';
        xCoord2plot5 = (xCoordMatAll5(indx2keep,dragTailStart:iFrame))';
        yCoord2plot5 = (yCoordMatAll5(indx2keep,dragTailStart:iFrame))';
        xCoord2plot6 = (xCoordMatAll6(indx2keep,dragTailStart:iFrame))';
        yCoord2plot6 = (yCoordMatAll6(indx2keep,dragTailStart:iFrame))';
        xCoord2plot7 = (xCoordMatAll7(indx2keep,dragTailStart:iFrame))';
        yCoord2plot7 = (yCoordMatAll7(indx2keep,dragTailStart:iFrame))';
        xCoord2plot8 = (xCoordMatAll8(indx2keep,dragTailStart:iFrame))';
        yCoord2plot8 = (yCoordMatAll8(indx2keep,dragTailStart:iFrame))';
        plotOrNot = 1;
    end
    
    %plot tracks
    if plotOrNot
        
        %plot basic tracks
        plot([xCoord2plot0(1,:); xCoord2plot0],[yCoord2plot0(1,:); yCoord2plot0],...
            'Color',[1 0.7 0.7],'LineWidth',4, 'Marker', 'o'); %light pink; the artificial repetition of the first line is for avoiding a mess in the first frame when tracks are not color-coded individually
                         
        %color individual tracks randomly if requested
        if colorTracks == 1
            plot(xCoord2plot1,yCoord2plot1,'Color',[1 0.7 0],'LineWidth',10); %orange
            plot(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',4, 'Marker', 'o'); %[1 0 0]
            plot(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',4, 'Marker', 'o'); %[0 1 0]
            plot(xCoord2plot4,yCoord2plot4,'Color','y','LineWidth',4, 'Marker', 'o'); %[1 1 0]
            plot(xCoord2plot5,yCoord2plot5,'Color','b','LineWidth',4, 'Marker', 'o'); %[0 0 1]
            plot(xCoord2plot6,yCoord2plot6,'Color','c','LineWidth',4, 'Marker', 'o'); %[0 1 1]
            plot(xCoord2plot7,yCoord2plot7,'Color','m','LineWidth',4, 'Marker', 'o'); %[1 0 1]
            plot(xCoord2plot8,yCoord2plot8,'Color',[0.6 0 1],'LineWidth',4, 'Marker', 'o'); %purple
        end
                
    end
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFramesMovie,movieVar,iFrame);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
end %(for iFrame = 1 : numFramesMovie)

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%% ~~~ end ~~~
