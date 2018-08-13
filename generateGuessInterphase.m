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
flagPlotAll = 1;
plot_overlay = 1;
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
imWienerM = imWiener .* maskForCell;

% run a steerable detector to locate line/curve features of a certain
% gaussian width
[res, theta, nms, ~] = steerableDetector( imWiener, 4,2);
res = res .* maskForCell;
nms = nms .* maskForCell;

% Basic Thresholding. Eliminate weak responses
T = multithresh( res, 1);
imRes = res;
imRes( res < T(1) ) = 0;
[a, ~,b] = steerableDetector( imRes, 4, 2);
T2 = multithresh(b, 2);
imB = b; imB( b < T2(1) ) = 0; imB = logical( imB);
imNms = bwareafilt( imB, [8, Inf] );

% define parameters for network reconstruction
costs.maxCost = 3;
costs.maxDistLink = 15;
costs.maxDistLinkPerp = 8;
costs.maxPhiDiff_EE = deg2rad( 60);
costs.maxPhiDiff_EL = deg2rad(60);
costs.EEvsELScalingFactor = 1.2;
costs.linkObjectSizeScalingFactor = 1.5;
costs.linkObjectSizeForScaling = 10;
params.costs = costs;
% Reconstruct MT network
imConn = reconstructNetwork(imNms, params);
imConn = bwareafilt( imConn, [15 Inf]);
% get individual lines/curves for reconstructed network
imHelper = findMinimalComponents(imConn, params);

if length( imHelper) > 5
    warning('Number of estimated lines is %d. Removing Small Lines', length( imHelper) )
    minLineSize = 20;
    rmLines = cellfun( @(v) sum( v(:) )<minLineSize, imHelper);
    imHelper( rmLines) = [];
end

fprintf('generateGuessInterphase: %d possible lines found.\n',2*length( imHelper)) 

%% Finding the MTOC (nucleation point of lines)

% It is also possible that the widest point along the line is an MTOC,
% and it has an overlap of tubes going different ways. We'll look for the
% brightest pixel in a convolved image.
imGauss = mat2gray( imgaussfilt( imageFull, 2) );
for jLine = 1 : length( imHelper)
    cLine = imHelper{jLine};
    cLine_dilated = imdilate( cLine, strel('disk', 1) );
    [~, MTOC(jLine).idx ] = max(imGauss(:) .* cLine_dilated(:) );    
    [MTOC(jLine).y,MTOC(jLine).x] = ind2sub( size(imageFull), MTOC(jLine).idx ); % note y and x are reversed because first index refers to row number, i.e the vertical axis
    
    % Lets also find the direction of nucleation for each MTOC
    st = regionprops( cLine, 'Orientation');
    MTOC(jLine).phi = deg2rad( st.Orientation) ;
end

% stop if no MTOCs found.
if length( imHelper) == 0
    MicrotubuleBank = [];
    disp('No MTOCs were found! No microtubules can exist for this cell.')
    return
end

% plot MTOCs on top of masked image
if flagPlotImp
    dispImg( imCell); hold on
    plot( [MTOC.x], [MTOC.y], 'b*', 'MarkerSize', 12, 'LineWidth', 3); hold off;
    title('MTOC locations')
    set(gcf, 'Name', 'mtoc_locations', 'NumberTitle', 'off')
end

% Now that we have the centroid locations, we can proceed with finding the
% lines

%% Estimating the lines

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

% For each Nucleation point, allow two microtubules oriented in opposite
% directions
stepSize = 2; visibility = 8; fieldOfVision = 60; %degrees

% Create edited image for estimate overlay
imDisp = mat2gray(imCell);
T = multithresh(imDisp( imDisp > 0), 3);
imDisp( imDisp < T(1)) = 0; 
imDisp = imadjust( mat2gray(imDisp), [T(1), 1]);

% Two microtubules per nucleation site. Assign helper image to each
% microtubule.
clear MicrotubuleBank
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
        helperImg =  mat2gray( imgaussfilt( imdilate( imHelper{ jNuc}, strel('disk', 1)) .* imCell) );
        MicrotubuleBank( jTube).helperImage = helperImg;
        MicrotubuleBank( jTube).display.image = imDisp;
        
    end
end

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
    plotEstimatedPointsAll( MicrotubuleBank( jTube) )
end

% Discard any dead microtubules
idxDead = find( [MicrotubuleBank.dead]);
MicrotubuleBank( idxDead) = [];
numTubes = length(MicrotubuleBank);
% Make Microtubule Bank into a class?

% Overlay estimated plot
if plot_overlay
    figName = sprintf( 'mt_estimates');
    figTitle = sprintf( 'Microtubule estimates');

    figure('Name', figName, 'NumberTitle', 'off')
    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
    subplot(121); imagesc( imDisp); axis equal; colormap gray; title('Original Image')
    set(gca, 'xlim', [1 size( imDisp ,1)], 'ylim', [1 size( imDisp ,2)]); set(gca, 'FontSize', 15);
    subplot(122); imagesc( imDisp); axis equal; colormap gray;
    set(gca, 'xlim', [1 size( imDisp ,1)], 'ylim', [1 size( imDisp ,2)]); hold on;
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

% plot all MT curves

% Overlay estimated plot
if plot_overlay
    figName = sprintf( 'mt_estimated_curves_%d', polyOrder);
    figTitle = sprintf( 'Microtubule estimated curves-%d', polyOrder);

    figure('Name', figName, 'NumberTitle', 'off')
    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
    subplot(121); imagesc( imDisp); axis equal; colormap gray; title('Original Image')
    set(gca, 'xlim', [1 size( imDisp ,1)], 'ylim', [1 size( imDisp ,2)]); set(gca, 'FontSize', 15);
    subplot(122); imagesc( imDisp); axis equal; colormap gray;
    set(gca, 'xlim', [1 size( imDisp ,1)], 'ylim', [1 size( imDisp ,2)]); hold on;
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


    function [imSeeds, nms, rotations] = findMTSeeds( imageOld)
        
        % Use a steerable filter to detect locations of lines (ridge detector)
        [res, theta, nms, rotations] = steerableDetector( imageOld,4,2.5, 180);

        % The suppressed response max has some noise in the background. We will use
        % the fact the most of the interphase lines are oriented along the
        % long-axis of the cell. We'll allow a wiggle room of 50 degrees on either
        % side. Then, we'll turn off pixels with lines oriented in disallowed
        % directions.
        % First, find the orientation of the cell.
        stats = regionprops( logicMask, 'Orientation', 'Centroid');
        phi = sign(stats.Orientation) * (90 - abs(stats.Orientation) );
        % Define line-angles that are allowed
        angRange = 50;
        badPixels = theta < deg2rad( phi - angRange) | theta > deg2rad( phi + angRange);

        % Turn off the bad pixels in the response map.
        imNMS = mat2gray(nms).*logicMask;
        imNMS( badPixels) = 0;

        % Now, our goal is two-fold: 1) find the number of lines, and 2) find the
        % brightest pixel on each line (which will act as the MTOC)

        % To achieve the first goal, we must distinguish between what is and what
        % isn't a line. We'll find a threshold.
        T = multithresh( imNMS(imNMS>0), 1);
        imThresh = imNMS; imThresh( imThresh < T) = 0;
        
        
        % Finally, as a fine-sieve, we'll use an area filter to remove anything
        % smaller than 5 pixels (too short of a line)
        imThresh = imThresh.* bwareafilt( logical(imThresh), [3 Inf]);

        imSeeds = imThresh;
        
    end

    function imTrees = growMTSeeds( imSeeds, rotations)
        
        % First, find the orientation of the cell.
        stats = regionprops( logicMask, 'Orientation', 'Centroid');
        phi = sign(stats.Orientation) * (90 - abs(stats.Orientation) );
        
        % find the closest index to the current cell orientation
        nAngs = size(rotations, 3);
        angs = [ linspace( 0, 180, nAngs/2), linspace( -180, 0, nAngs/2) ];
        [~, idx] = min( abs(angs-phi) );
        indices = mod( idx: idx+nAngs, nAngs);

        indices = [indices, flip(indices) ];
        indices( indices == 0) = 1;

%         figure('WindowState', 'maximized')

        imComp = nms .* logicMask;
        T1 = multithresh( imComp( imComp > 0), 2);

        st1 = regionprops( logical(imSeeds), 'PixelIdxList', 'Orientation');
        phi0 = [st1.Orientation];
        thDif = 20; %degrees

        imAdd = logical( imSeeds);
        for j = 1 : length(indices)

            currIdx = indices(j);
%             currAng = angs( currIdx);

            imCurr = rotations(:,:,currIdx) .* logicMask;
            T2 = multithresh( imCurr(imCurr > 0), 2);
            imT = imCurr; imT( imT < T2(1) | imT < T1(1) ) = 0;

%             subplot(121); imagesc(imT); colormap gray; axis equal

            st2 = regionprops( logical(imT), 'PixelIdxList', 'Orientation');

            for jj = 1 : length(st2)

                overlap = 0; orient_match = 0;
                % is there any overlap of this object?
                imObj = 0*imSeeds;
                imObj( st2( jj).PixelIdxList ) = 1;
                if sum( imAdd(:) & imObj(:) ) > 1; overlap = 1; end

                % do the orientations match?
                oldOr = phi0;
                newOr = st2( jj).Orientation;
                if abs( abs(oldOr) - abs(newOr) ) < thDif; orient_match = 1; end
                if orient_match && overlap

                    imAdd = logical( imAdd + imObj);
                    phi0 = st2( jj).Orientation;
                else
                    st3 = regionprops( imAdd, 'Orientation');
                    phi0 = [st3.Orientation];
                end

            end


        end
        imTrees = imAdd;
        
    end


end

