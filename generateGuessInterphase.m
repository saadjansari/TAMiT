function MicrotubuleBank = generateGuessInterphase(imageForGuess, maskForCell, SavePath)
% generateGuessInterphase: Generates a guess for the microtubules in
% interphase.
%   For the first argument, it is possible to give in a 3D array where the
%   third dimension corresponds to time. In that case, the program will use
%   a Maximum Intensity Projection in time to obtain a 2D image. Before
%   using multiple times, please confirm that the imaging time-scale is
%   smaller than the dynamic time-scale.


% plot flags for visual control/debugging
flagPlotImp = 1;
flagPlotAll = 0;
if flagPlotAll
    flagPlotImp = 1;
end

% If image2D is a 3D array, take a MIP
if size( imageForGuess, 3) > 1
    warning(['generateGuessInterphase: The user provided a 3D array for the first argument. ',...
        '3rd dimension is assumed as time and an MIP is being taken'])
    imageFull = max(imageForGuess, [], 3);
end

% Apply the mask to the full image
imCell = imageFull .* maskForCell;

% Filter out noise using the wiener2 filter (evens out variations in local
% standard deviaitons)
imWiener = wiener2( imageFull, [3 3]);

% We need to get the skeleton of the microtubules.
% We will proceed in the following manner:
% 1. The original cell image is noisy, so we will apply a gaussian filter
%    to it equal to the PSF of our microscope of reduce the effects of shot
%    noise. We will alter the image to remove the background for the
%    estimation process.
% 2. Next, we use a steerable detector to find line-responses.
% 3. We say that a response > resLine corresponds to an actual line, and
%    other pixels do not contain lines. This is a thresholding step. We
%    find the threshold via otsu thresholding however the threshold value
%    is bounded from below.
% 4. We find the single-pixel-width lines that best explain the response
%    map.

imEstimate = imageFull; imEstimate( imEstimate == 0) = median( imEstimate( imEstimate ~= 0) );
T = multithresh(imEstimate( maskForCell), 3);
imEstimate( imEstimate < T(1)) = 0;
imEstimate = imadjust( mat2gray(imEstimate), [T(1), 1 ] );
imGaussian = mat2gray( imgaussfilt( imEstimate, 1.5) );
imGaussianMasked = mat2gray( imGaussian .* maskForCell);

[~, ~, nonmax] = steerableDetector( imGaussian, 4, 2); SWL = nonmax .* maskForCell;
T3 = multithresh( SWL( SWL > 0), 4); SWL( SWL < T3(1) ) = 0; SWL = logical(SWL);

imBombed = bwareafilt( SWL, [10, Inf]); imBombed = bombNetworkJunctions( imBombed);
imBombed = bwareafilt( imBombed, [10, Inf]);

% cause breaks in the skeleton whenever there is a sharp orientation change
% in an object.
% We need a function that can help us find local orientations of each
% object in the skeleton
skelProps = analyzeSkeleton( imBombed);
xRm = []; yRm = [];
for jj = 1 : length(skelProps)
    signal = smooth( skelProps(jj).thetaRaw, 3 );
    
    s1 = signal; s2 = signal;
    s1(end-1:end) = []; s2( 1:2)=[];
    x = skelProps(jj).orderedPixelList(:, 2); x(1) = []; x(end) = [];
    y = skelProps(jj).orderedPixelList(:, 1); y(1) = []; y(end) = [];
    deriv3 = s1-s2;
    idxKill = find( abs(deriv3) > 0.3);
    xRm = [xRm; x(idxKill) ];
    yRm = [yRm; y(idxKill) ];
end

idxKill = sub2ind( size(imBombed), yRm, xRm);
imTemp = imBombed;
imTemp( idxKill) = 0;
imBroken = bwareafilt( imTemp, [5, Inf]);

stD = regionprops( imdilate(imBroken, strel('disk', 1) ), imGaussian, 'MeanIntensity', 'PixelIdxList');
st = regionprops( imBroken, imGaussian, 'MeanIntensity', 'PixelIdxList');
idxAvoid = [ cat( 1, stD(:).PixelIdxList); find(maskForCell == 0) ];
vals = imGaussian(:); vals( idxAvoid) = NaN;
mu = median( vals( ~isnan(vals) ) );
sig = std( vals( ~isnan(vals) ) );
objKeep = find([st.MeanIntensity] > mu+sig);
imNew = 0*imCell;
for jj = 1 : length(objKeep)
    imNew( st( objKeep(jj) ).PixelIdxList ) = 1;
end

im2 = repmat( mat2gray( imGaussian .* maskForCell).* ~SWL, 1, 1, 3); im2(:,:,3) = 0*imCell; im2( :,:,2) = SWL;
im3 = repmat( mat2gray( imGaussian .* maskForCell).* ~imNew, 1, 1, 3); im3(:,:,3) = 0*imCell; im3( :,:,2) = imNew;


% Now lets reconstruct the network by connecting any obviously connected
% pieces. We put a strick constraint on orientation matching but allow
% linking to objects as far as 15 pixels (~1.5 micron)

% define parameters for network reconstruction
costs.maxCost = 3;
costs.maxDistLink = 20;
costs.maxDistLinkPerp = 2;
costs.maxPhiDiff_EE = deg2rad( 20);
costs.maxPhiDiff_EL = deg2rad(10);
costs.EEvsELScalingFactor = 1.5;
costs.linkObjectSizeScalingFactor = 1.2;
costs.linkObjectSizeForScaling = 10;
paramsReconstruct.costs = costs;

clear costs
costs.maxDistLink = 12;
costs.maxPhiDiff_EE = deg2rad( 70);
paramsFindComp.costs = costs;

% Reconstruct MT network
imConn = reconstructNetwork(imNew, imCell, paramsReconstruct);

im4 = repmat( mat2gray( imGaussian .* maskForCell).* ~imConn, 1, 1, 3); im4(:,:,3) = 0*imCell; im4( :,:,2) = imConn;
dispImg( imCell, im2, im3, im4, [ 2 2]); set(gcf, 'Name','cleaned_up_network', 'NumberTitle', 'off')
imHelper = findMinimalComponents( logical(imConn), paramsFindComp);

rmSmallLines = 1;

if rmSmallLines
    fprintf('Number of estimated lines is %d. Removing Small Lines\n', length( imHelper) );
    minLineSize = 10;
    rmLines = cellfun( @(v) sum( v(:) )<minLineSize, imHelper);
    imHelper( rmLines) = [];
end
fprintf('generateGuessInterphase: %d possible tracks found.\n',length( imHelper)) 

%% Finding the MTOC (nucleation point of microtubules)

% stop if no MTOCs found.
if length( imHelper) == 0
    MicrotubuleBank = [];
    disp('No Tracks and no MTOCs were found! No microtubules can exist for this cell.')
    return
end

% It is also possible that the widest point along the line is an MTOC,
% and it has an overlap of tubes going different ways. We'll look for the
% brightest pixel in a convolved image.

imGauss = mat2gray( imgaussfilt( imageFull, 2) );
for jLine = 1 : length( imHelper)
    
    cLine = imHelper{jLine};
%     cLine_dilated = imdilate( cLine, strel('disk', 1) );
    [~, MTOC(jLine).idx ] = max(imGauss(:) .* cLine(:) );  
    [MTOC(jLine).y,MTOC(jLine).x] = ind2sub( size(imageFull), MTOC(jLine).idx ); % note y and x are reversed because first index refers to row number, i.e the vertical axis
    
    % Lets also find the direction of nucleation
    boneProps = analyzeBone( cLine);
    MTOC(jLine).phi = boneProps.thetaRaw(boneProps.orderedIdxList == MTOC(jLine).idx);
    
end

% plot MTOCs on top of masked image
if flagPlotImp
    imRGB = mat2gray(repmat( imGaussianMasked, 1, 1, 3)); imRGB( :,:, 2:3) = 0;
    dispImg( imRGB); hold on
    plot( [MTOC.x], [MTOC.y], 'w*', 'MarkerSize', 12, 'LineWidth', 3); hold off;
    title('MTOC locations')
    set(gcf, 'Name', 'mtoc_locations', 'NumberTitle', 'off')
end

%% Estimating the lines

clear MicrotubuleBank
% To estimate the lines, we'll take the following approach:
% 1) Start at the MTOC (allow 2 nucleations on either side)
% 2) For each nucleation, we'll find the angle that gives the maximum in a
%    radially integrated signal, up to some max radius( which we will call
%    the visibility)
% 3) We'll propagate from the MTOC to a new point along the optimum angle,
%    which a defined distance away from the MTOC (this we will call the
%    step size). Throughout this process, we might look for optimum angles
%    within a certain range related to the orientation of the Tube (this we
%    will refer to as the field of vision).
% 4) We will continue propagating until we can't find an optimum angle
%    anymore (all our surroudings appear the same. This will be based on
%    some optimality condition. At this point we will change our step size
%    to half its original value and see if we can propogate a smaller
%    distance forward.
% 5) This will conclude our initial estimation phase.
% 6) We will fit a polynomial through each initial estimate (quadratic or
%    cubic). This will be our final estimate for the microtubules.

% Ensure 1 MTOC, 2 microtubules per track. Assign helper images and correct
% phi to each of the 2 microtubules per track.
for jNuc = 1 : length(MTOC)
    imt = 1 + 2*(jNuc-1);
    fmt = imt+1;
    for jTube = imt:fmt
        if mod( jTube, 2)
            MicrotubuleBank( jTube) = Microtubule( imCell, MTOC( jNuc).x, MTOC( jNuc).y, MTOC( jNuc).phi );
        else
            MicrotubuleBank( jTube) = Microtubule( imCell, MTOC( jNuc).x, MTOC( jNuc).y, MTOC( jNuc).phi - pi );
        end
        % Assign the helper image
        helperImg =  mat2gray( imgaussfilt( 1.0*imdilate( imHelper{ jNuc}, strel('disk', 1)) .* imCell ) );
        MicrotubuleBank( jTube).helperImage = helperImg;
        MicrotubuleBank( jTube).display.image = imGaussianMasked;
        
    end
end


% For each Nucleation point, allow two microtubules oriented in opposite
% directions
stepSize = 2; visibility = 8; fieldOfVision = 60; %degrees

% Assign identification and color to each microtubule. 
numTubes = length( MicrotubuleBank);
colors = distinguishable_colors( numTubes, {'w', 'k'} );
for jTube = 1 : numTubes
    MicrotubuleBank( jTube).id = jTube;
    MicrotubuleBank( jTube).display.color = colors( jTube, :);
end

% Estimate MT locations
for jTube = 1 : numTubes
    MicrotubuleBank( jTube) = EstimateMicrotubulePoints(  MicrotubuleBank( jTube), stepSize, visibility, fieldOfVision);
%     plotEstimatedPointsAll( MicrotubuleBank( jTube) )
end


% Discard any dead microtubules
idxDead = find( [MicrotubuleBank.dead]);
MicrotubuleBank( idxDead) = [];
numTubes = length(MicrotubuleBank);

% Make Microtubule Bank into a class?

% Overlay estimated plot
if flagPlotAll
    figName = sprintf( 'mt_estimates');
    figTitle = sprintf( 'Microtubule estimates');

    figure('Name', figName, 'NumberTitle', 'off')
    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
    subplot(121); imagesc( imGaussianMasked); axis equal; colormap gray; title('Original Image')
    set(gca, 'xlim', [1 size( imGaussianMasked ,1)], 'ylim', [1 size( imGaussianMasked ,2)]); set(gca, 'FontSize', 15);
    subplot(122); imagesc( imGaussianMasked); axis equal; colormap gray;
    set(gca, 'xlim', [1 size( imGaussianMasked ,1)], 'ylim', [1 size( imGaussianMasked ,2)]); hold on;
    LH = []; L = {};
    for jTube = 1 : numTubes

        x0 = MicrotubuleBank( jTube).estimatedPoints(1,:);
        y0 = MicrotubuleBank( jTube).estimatedPoints(2,:);
        c = MicrotubuleBank( jTube).display.color;
        lw = MicrotubuleBank( jTube).display.LineWidth+2;
        ms = MicrotubuleBank( jTube).display.MarkerSize;
        id = MicrotubuleBank( jTube).id;

        % plot all estimated points with links in between them
        plot( x0, y0, '-', 'Color', c, 'LineWidth', lw, ...
            'Marker', '.', 'MarkerSize', ms)

        % plot mtoc location
        plot( x0(1), y0(1), 'r*', 'LineWidth', lw/2, ...
            'MarkerSize', ms+4)

        % create legend entry
        LH(jTube) = plot(nan, nan, '-', 'Color', c, 'MarkerSize', ms, 'LineWidth', lw); L{jTube} = ['MT ', num2str(id)];    

    end
    hold off
    legend(LH, L);
    set(gca, 'FontSize', 15);
    title( figTitle);
end

% Estimate polynomial curve coefficients
polyOrder = 3;
for jTube = 1 : numTubes
    MicrotubuleBank( jTube) = EstimateMicrotubuleCurve( MicrotubuleBank( jTube), polyOrder);
    PlotMicrotubuleCurve( MicrotubuleBank( jTube), 'estimate')
end

% Overlay estimated curve plot
if flagPlotImp
    figName = sprintf( 'mt_estimated_curves_%d', polyOrder);
    figTitle = sprintf( 'Microtubule estimated curves-%d', polyOrder);

    figure('Name', figName, 'NumberTitle', 'off')
    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
    subplot(121); imagesc( imGaussianMasked); axis equal; colormap gray; title('Original Image')
    set(gca, 'xlim', [1 size( imGaussianMasked ,1)], 'ylim', [1 size( imGaussianMasked ,2)]); set(gca, 'FontSize', 15);
    subplot(122); imagesc( imGaussianMasked); axis equal; colormap gray;
    set(gca, 'xlim', [1 size( imGaussianMasked ,1)], 'ylim', [1 size( imGaussianMasked ,2)]); hold on;
    LH = []; L = {};
    for jTube = 1 : numTubes

        
        t = linspace(0,1);
        px = MicrotubuleBank( jTube).estimatedCoef(1,:);
        py = MicrotubuleBank( jTube).estimatedCoef(2,:);
                
        x = polyval( px, t );
        y = polyval( py, t );
        c = MicrotubuleBank( jTube).display.color;
        lw = MicrotubuleBank( jTube).display.LineWidth+2;
        ms = MicrotubuleBank( jTube).display.MarkerSize;
        id = MicrotubuleBank( jTube).id;

        % plot all estimated points with links in between them
        plot( x, y, '-', 'Color', c, 'LineWidth', lw, ...
            'Marker', '.', 'MarkerSize', ms)

        % plot mtoc location
        plot( x(1), y(1), 'r*', 'LineWidth', lw/2, ...
            'MarkerSize', ms+4)

        % create legend entry
        LH(jTube) = plot(nan, nan, '-', 'Color', c, 'MarkerSize', ms, 'LineWidth', lw); L{jTube} = ['MT ', num2str(id)];    

    end
    hold off
    legend(LH, L);
    set(gca, 'FontSize', 15);
    title( figTitle);
end

% Save all figures
SaveOpenFigsToFolder( SavePath, 'png', 1)

    function networkout = bombNetworkJunctions( networkIn)
        
        
        % bomb places where their are any connections
        
        % clear the border
        networkIn( 1, :) = 0; networkIn( :, 1) = 0; networkIn( end, :) = 0; networkIn( :, end) = 0;
        networkout = networkIn;
        cc = bwconncomp( networkout); nObj = cc.NumObjects;
        countBombed1 = 0; countBombed2 = 0; countBombed3=0;
        for jObj = 1 : nObj

            % get all its pixels, and turn on just the object
            pixList = cc.PixelIdxList{ jObj};
            imObj = 0 * networkout;
            imObj( pixList) = 1;
            % get cartesian indices of pixel
            [yList, xList] = ind2sub( size(imObj), pixList );

            % check if any pixel satisfies our criterion for being a connection
            % pixel
            for jPix = 1 : length(pixList)
                if pixList( jPix) == 11400
                    stopH = 1;
                end
                y = yList( jPix); x = xList( jPix);

                % get the local 3x3 neighborhood
                nhood = imObj( y-1:y+1, x-1:x+1);

                nhood( 5) = 0;
                % apply criterion for endpoint
                if sum( nhood( :)) == 3

                    idx = find( nhood);
                    idx = idx( idx ~= 5); % not the current pixel
                    [yy, xx] = ind2sub( size(nhood), idx);
                    dist = [];
                    for jIdx = 1 : length(idx)
                        for kIdx = 1: length(idx)
                            if jIdx ~= kIdx
                                dist = [ dist, norm( [yy(jIdx)-yy(kIdx), xx(jIdx)-xx(kIdx)] ) ];
                            end
                        end
                    end
                    if all( dist >= 1.4)
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
                        countBombed1 = countBombed1+1;
                    end

                elseif sum( nhood(:)) == 4 && sum( nhood(2:2:end) ) ~= 4

                    idx = find( nhood);
                    idx = idx( idx ~= 5); % not the current pixel
                    [yy, xx] = ind2sub( size(nhood), idx);
                    Mindist = 10*ones( 1, length(idx) );
                    for jIdx = 1 : length( idx)
                        for kIdx = 1: length( idx)
                            if jIdx ~= kIdx
                                dist = norm( [yy(jIdx)-yy(kIdx), xx(jIdx)-xx(kIdx)] );
                                Mindist( jIdx) = min( [ Mindist( jIdx), dist ] );
                            end
                        end
                    end
                    if sum( Mindist >= 2) >= 2 || ( sum( Mindist >= 2) >= 1 && sum( Mindist >= 1.4) >= 1)
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
                        countBombed2 = countBombed2+1;
                    end
                elseif sum( nhood(:))  == 5
                    idx = find( nhood);
                    idx = idx( idx ~= 5); % not the current pixel
                    [yy, xx] = ind2sub( size(nhood), idx);
                    Mindist = 10*ones( 1, length(idx) );
                    for jIdx = 1 : length( idx)
                        for kIdx = 1: length( idx)
                            if jIdx ~= kIdx
                                dist = norm( [yy(jIdx)-yy(kIdx), xx(jIdx)-xx(kIdx)] );
                                Mindist( jIdx) = min( [ Mindist( jIdx), dist ] );
                            end
                        end
                    end
                    if ( ( sum(Mindist == 1) == 5) && ( ~istriu(nhood) ) && ( ~istril(nhood) ) ) ||...
                            ( ( sum( Mindist == 1) == 4) && ( sum( Mindist >= 1.3) == 1) ) ||...
                            ( sum( Mindist >= 2) == 2)
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
                        countBombed3 = countBombed3+1;
                    end
                elseif sum( nhood(:))  >= 6
                    [y,x] = ind2sub( size(networkout), pixList( jPix) );
                    networkout( y-1:y+1, x-1:x+1 ) = 0;
                    countBombed3 = countBombed3+1;
                end        

            end

        end
        
    end

    function boneProps = analyzeBone( bone)
        % findCurveVariance: find the variance of the local curvatures of
        % an object as well as some other properties

        % find endpoint of line
        bone = logical( bone);
        bw = bwmorph( bone, 'endpoints');
        endPtIdx = find(bw == 1);
        [ y1, x1] = ind2sub( size(bone), endPtIdx);
        
        % pick an endpoint and measure geodesic distance to each other
        % pixel to get order list of connected pixels
        mat_dist = bwdistgeodesic( bone, x1(1), y1(1), 'quasi-euclidean'); %'quasi-euclidean' for 8-connectivity
        comp = find( bone);
        comp(:,2) = mat_dist(comp(:,1));
        ordered_list_ind = sortrows(comp,2);
        [ordered_list_sub(:,1) ordered_list_sub(:,2)] = ind2sub( size(bone), ordered_list_ind(:,1) ); 
        
        xyList = ordered_list_sub;
        numPix = size( xyList,1);
        % we'll use a moving filter of 5 pixels to find the local
        % orientations
        filtSize = 5;
        filtHalf = (filtSize-1)/2 * mod(filtSize, 2) + filtSize/2 * ~ mod(filtSize, 2);

        thetaList = [];
        for jPix = 1 : numPix
            
            if jPix <= filtHalf
                startPix = 1;
            else
                startPix = jPix - filtHalf;
            end
            if jPix > numPix - filtHalf
                endPix = numPix;
            else
                endPix = jPix + filtHalf;
            end
            
            localLine = xyList( startPix : endPix, :)';
            
            yD = localLine(1, end) - localLine(1, 1);
            xD = localLine(2, end) - localLine(2, 1);
            
%             xD = diff(localLine(1,:) );
%             yD = diff(localLine(2,:) );
            
            thetaNextPossible = mean( atan2(yD, xD) ) + [-2*pi, 0, +2*pi];
            if ~isempty( thetaList)
                [~, thetaNextIdx] = min( abs( thetaNextPossible - thetaList(end) ) );
            else
                thetaNextIdx = 2;
            end
            thetaList = [ thetaList, thetaNextPossible( thetaNextIdx) ];
            
        end
        
        % store information about line/curve
        boneProps.length = numPix;
        boneProps.theta = mean( thetaList);
        boneProps.thetaVar = std( thetaList);
        boneProps.thetaRaw = thetaList;
        boneProps.orderedPixelList = xyList;
        boneProps.orderedIdxList = ordered_list_ind(:,1);
    end

    function skeletonProps = analyzeSkeleton( skeleton)
        % analyzeLineImages: runs the findLineImageVariance function on
        % each image inside the imageCellArray and stores in the lineProps
        % structure
        
        cc = bwconncomp( skeleton);
        for jBone = 1 : cc.NumObjects
            
            imBone = 0*skeleton;
            imBone( cc.PixelIdxList{jBone} ) = 1;
            
            skeletonProps(jBone) = analyzeBone( imBone);
            
        end
        
    end

end

