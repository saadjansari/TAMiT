function [coords, orients] = searchForHighIntensityLines3( imGauss, imLines, imMask, estParams)

    % To estimate the lines, we'll take the following approach:
    % 1) Start at the brightest point in the image, and find the optimal orientation of any underlying lines  by radially integrating intensity.
    % 2) Once we have a start point and a direction, we continue searching for the underlying line by looking for intensity peaks until a threshold is met.
    % 3) We take this estimated microtubule, and the pixel region it spans and set it to the median value of the original image. We say that we have removed the intensity due to this microtubule.
    % 4) We repeat steps 1-3 until the microtubule estimates are smaller than a certain length
    % 5) We will fit a polynomial through each initial estimate (quadratic or
    %    cubic). This will be our final estimate for the microtubules.


    % plotflags {{{
    plotflag_iter = 0;
    plotflag_zslice = 0;
    plotflag_costcalc = 0;
    plotflag_link = 0;
    plotflag_allZ = 0;
    dispflag_imageProfile = 0;
    plotflag_finalLinkedFilaments = 0;
    % }}}

    % Extract Parameters {{{
    % extract parameters
    vis = estParams.visibility;
    fov = estParams.fieldOfView;
    step = estParams.stepSize;
    minMTlength = estParams.minMTlength;
    maxIter = estParams.maxMTsearchIterations;
    costAccept = estParams.cost_accept_link;
    maxFails = 2; 
    imEst = imLines.*imMask;
    % }}}

    % Parse Image to Parameters {{{
    % Get pixel intensity distribution information from the intensity image
    imMask3 = logical( 0*imGauss);
    for jj=1:size(imGauss,3), imMask3(:,:,jj) = logical( imMask); end

    imVals = imGauss( imMask3);
    med = median( imVals );
    otsu = multithresh( imVals, 2); otsu = otsu(2);
    mu = mean( imVals);
    sig = std( imVals);
    muMin = min( [mu, med]);
    
    if dispflag_imageProfile
        disp( sprintf( 'Mean = %.2f', mu) );
        disp( sprintf( 'STD = %.2f', sig) );
        disp( sprintf( 'Median = %.2f', med) );
        disp( sprintf( 'Otsu = %.2f', otsu) );
    end
    % }}}

    % Z-loop parameters {{{
    % initialize parameters for looping
    intThreshMTOC = mu+3*sig;
    intThreshMT = mu+sig;
    disp( sprintf( ' MTOC Intensity Threshold = %.2f', intThreshMTOC) )
    disp( sprintf( ' MT Intensity Threshold = %.2f', intThreshMT) )
    intMax = 1.2*intThreshMTOC;
    lineIdx = 1; 
    clabel = 1;
    % }}}

    % We will loop over all frames in Z and make detections of microtubules
    for zFrame = 1 : size(imGauss, 3)
        
        % Init Z-specific params and while-loop params {{{
        disp( sprintf('  Z-Slice = %d', zFrame ) )
        imEst = imLines(:,:,zFrame).*imMask3(:,:,zFrame);
        imEstG = imEst.*imGauss(:,:,zFrame);
        imGaussZ = imGauss(:,:,zFrame);
        fails = 0;
        iter = 0;
        coords = {};
        orients = [];
        currIndex = 0;
        LineNum = [];
        mtocZ = [];
        mtocZint = [];

        intThreshMTOC = mean( imGaussZ( find(imMask) ) ) + 2*std( imGaussZ( find(imMask) ) ); 
        intThreshMT = mean( imGaussZ( find(imMask) ) );
        intThreshMTRatio = 0.7;
        disp( sprintf( 'intMTOC_Z = %.2f', intThreshMTOC) )
        disp( sprintf( 'intMT_Z = %.2f', intThreshMT) )
        % }}}

        % Find microtubules in a single 2D frame {{{
        while (intMax > intThreshMTOC || max( logical( imEstG(:) ) ) > intThreshMTOC) && iter < maxIter && fails < maxFails
             
            iter = iter +1;
            
            % Find Bright Point in Image
            [ intMax, idxMax] = max( imEst(:) );
            [ yMax, xMax] = ind2sub( size(imEst), idxMax);

            % Does MT pass through this Bright Point? 
            [coordTemp, success, orientStart] = estimateNextPoint( [ xMax; yMax], 0, step, vis, 540, imEst);
            
            
            % if MT passes, search for MT on both sides 
            if success
                
                currIndex = currIndex+1;
                % 2 angles because 2 sides
                angOrient = [ orientStart, orientStart+pi]; 
                for jAng = 1:2
                    ang = angOrient( jAng);
                    coordMT{jAng} = estimateMicrotubulePoints( [xMax; yMax], ang, imEst, step, vis, fov, 1);
                end

                coordFull = [ flip( coordMT{1}, 2), coordMT{2} ]; coordFull ( :, size(coordMT{1}, 2) ) = [];

                % find real MTOC using the intensity image informaiton instead of the steerable filtered image
                [cX,cY] = interpolateCoords( coordFull(1,:), coordFull(2,:), round(4*step) ); 
                indXY = sub2ind( size( imGaussZ), round(cY), round(cX) );
                [intMax, mtocIdx] = max( imGaussZ( indXY) );
                [ mtoc(2), mtoc(1)] = ind2sub( size(imGaussZ), indXY(mtocIdx) );
                disp( sprintf( '   iter = %d, int_MTOC = %.2f, int_MT = %.2f' , iter, intMax, mean(imGaussZ(indXY)) ) )
                mtocZ = [mtocZ, [mtoc(1);mtoc(2)]];
                mtocZint = [ mtocZint, intMax];

                % append found MTs if min length bigger than some threshold and intensityMTOC is biggger than threshold
                if size(coordFull, 2)*step > 1.5*minMTlength && intMax > intThreshMTOC && ( sum(imGaussZ( indXY) > intThreshMT)/ length(indXY) ) > intThreshMTRatio
                    addedMT = 1;
                    coords = {coords{:}, coordFull};
                else, addedMT = 0; end

%                 if addedMT && sum( imGaussZ(indXY) > intThreshMT)/length( indXY) > 0.7
%                     addedMT = 1; coords = {coords{:}, coordFull};
%                 end

                % Remove the MT from the image, and de-gauss the steetable MTOC
                imPost = removeMicrotubuleFromImage( imEst, coordFull, [xMax, yMax], 3, 0, 0, imMask);
                
                % display iteration plot if flag activated
                displayIterationPlot( imGaussZ, imEst, imPost, imMask, coordFull, [xMax, yMax], mtoc, iter, plotflag_iter); 
                if plotflag_iter, set(gcf, 'NumberTitle', 'off', 'Name', sprintf('z=%d_iter=%d' , zFrame, iter) ), end
                if addedMT && ~isempty(LineNum), LineNum( end+1) = LineNum(end)+1; elseif addedMT && isempty(LineNum), LineNum(1) = addedMT; end 

            else
                % de-gauss the steetable MTOC
                imPost = removeMicrotubuleFromImage( imEst, [], [xMax, yMax], 3, 0, 0, imMask);
                fails = fails+1;
                disp( sprintf( '   iter = %d , fails = %d', iter, fails) )
            end
            imEst = imPost;
        end
        
        % }}}
     
        % Save information about detected microtubules for use in connecting MTS from different z-slices {{{
        % these microtubules are objects that should be assigned a unique label, coordinates, an MT axis, a z-frame number, and some intensity
        for jC = 1 : length(coords)
            mts( clabel).coords = coords{jC};
            mts( clabel).z = zFrame;
            mts( clabel).int = mtocZint(jC);
            mts( clabel).int = mtocZ( :, jC);
            clabel = clabel+1;
        end
        % }}}

        % Z_plot {{{
        lineIdx = displayZplot( imGaussZ, imLines, imMask, zFrame, coords, LineNum, lineIdx, plotflag_zslice);
        if plotflag_zslice, set( gcf, 'NumberTitle', 'off', 'Name', sprintf('z=%d', zFrame) ), end
        % }}}

    end
    

    % All Z-plot {{{
    displayAllZplot( mean(imGauss, 3), mts, plotflag_allZ);
    if plotflag_allZ && exist( 'obj', 'var')
        set(gcf, 'Name', sprintf('cell_%d_allZ' , obj.sourceInfo.cellNumber) )
    end
    % }}}
        
    % Analyze each detected line and save its information
    mts = analyzeLine( mts, imMask);

    % find link costs between microtubules in structure
    [costs, mts, linkArray] = findLinkCosts( mts, costAccept, plotflag_costcalc); 
    % link mts if the cost is below a certain threshold
    [imagesMT, imagesID, mts] = linkMts( mts, costs, costAccept, imMask, linkArray, plotflag_link);
    
    [coords, orients] = findFinalMT( imagesMT, imagesID, mts, imGauss, minMTlength, plotflag_finalLinkedFilaments);

    if isfield(estParams, 'figureFlag') && estParams.figureFlag == 1
        sP = [pwd, filesep, 'PaperFigures/Estimation/'];
        save( [sP, 'mts&imgs.mat'], 'mts', 'imagesMT')
    end
    
    % findRealMTOCViaIntensityImage {{{
    function [coordMT, orients, intMTOC, mtoc] = findRealMTOCViaIntensityImage( coordMT, imInt);

        % get a single conencted line from the two detected microtubules. Use the true intensity information along the detected line to determine an ideal location for the MTOC.
        coordsAll = round( [flip( coordMT{1}, 2), coordMT{2}(:, 2:end)]);
        indAll = sub2ind( size(imInt), coordsAll(2,:), coordsAll(1,:) );

        [intMTOC, indMax] = max( imInt( indAll) );
        mtoc = coordsAll( :, indMax);
        coordMT{1} = coordsAll(:, indMax:-1:1);
        coordMT{2} = coordsAll(:, indMax:1:end);
       
        for jmt = 1 : 2
            try
            delX = diff( coordMT{jmt}( 1, 1:2));
            delY = diff( coordMT{jmt}( 2, 1:2));
            catch, delX=NaN; delY=NaN; end
            orients(jmt) = atan2( delY,delX);
        end

    end
    %  }}}

    % appendMTifConditionsMet {{{
    function [coords, orients, addedMT] = appendMTifConditionsMet( coords, coordMT, orients, orientMT, minMTlength, step, intMax, intThresh, lineNumber)
        
        addedMT = [];
        for jAng = 1 : 2
            if size( coordMT{jAng}, 2)*step > minMTlength && intMax >= intThresh
                addedMT = [addedMT, lineNumber];
                coords = { coords{:}, coordMT{jAng} };
                orients = [orients, orientMT(jAng)];
            end
        end

    end
    %  }}}
    
    % removeMicrotubuleFromImage {{{
    function imPost = removeMicrotubuleFromImage( imPrior, coordMT, deGaussPoint, deGaussPointSpread, mtVal, ptVal, imMask)
        
        imPost = imPrior;
        imPost( imPost < 0) = 0;

        if isempty( coordMT), skipMTremoval = 1; else, skipMTremoval = 0; end

        if ~skipMTremoval
            imMT = logical( 0*imPrior);
            [cdX, cdY] = interpolateCoords( coordMT(1,:), coordMT(2,:), 10);
            idxMT = sub2ind( size(imPrior), round(cdY), round( cdX) );
            imMT( idxMT) = 1;
            imMT = imdilate( imMT, strel('disk', 2) );
            idxMT = find( imMT);
            imPost( idxMT) = 0;
            minV = min( imPost(:) ); maxV = max( imPost(:) ); 
            imPost = mat2gray( imPost) * ( maxV-minV) + minV;
            imPost = imPost .* imMask;
            imPostLog = logical(imPost); 
            imPost = bwareafilt( imPostLog, [15, Inf]) .* imEst;
        end

        % de-intensity the mtoc point
        for jX = 1 : size(imPost, 2); for jY = 1 : size(imPost, 1)
            if norm( [jX, jY] - deGaussPoint ) <= deGaussPointSpread
                imPost(jY, jX) = 0;
            end
        end; end
        imPost = imPost.*imMask;

    end
    %  }}}

    % displayIterationPlot {{{
    function displayIterationPlot( imOriginal, imPrior, imPost, imMask, coordMT, mtocFake, mtocReal, iter, plotflag); 

        if ~plotflag, return, end

        figure('NumberTitle', 'off', 'Name', ['Iter_' num2str(iter)] )
        subplot(221)
        imagesc(imOriginal.*imMask); colormap gray
        title('Original Image')
        subplot(222)
        imagesc( imPrior); colormap gray; hold on; plot( mtocFake(1), mtocFake(2), 'b*', 'MarkerSize', 14); hold off;
        title( 'Search Image: MTOC' )
        subplot(223)
        imagesc( imPrior); colormap gray; hold on; plot( coordMT(1, :), coordMT(2,:), 'r-', 'LineWidth', 2); plot( mtocReal(1), mtocReal(2), 'b*', 'MarkerSize', 14); hold off
        title( 'Microtubule found' )
        subplot(224)
        imagesc(imPost); colormap gray; title('Image for next iteration')

    end
    % }}}

    % displayZplot {{{
    function nextIdx = displayZplot( im2D, imLines, imMask, currZframe, coords, lineNum, startIdx, plotflag)
        
        if ~plotflag, nextIdx = startIdx+length(coords); return, end
        
        LH=[]; L={};
        fS = size(im2D, 1);
        figure;
        subaxis(1,3,1); imagesc( im2D.*imMask ); colormap gray; axis equal; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
        subaxis(1,3,2); imagesc( imLines(:,:,currZframe).*imMask ); colormap gray; axis equal; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
        subaxis(1,3,3); imagesc( im2D.*imMask ); colormap gray; axis equal; hold on; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
        
        cNum = lineNum;
        lineNum = lineNum+startIdx-1;

        col = distinguishable_colors( max(cNum), {'k', 'w'});
        for jmt = 1 : length(coords)
            plot( coords{jmt}(1,:), coords{jmt}(2,:), 'color', col( cNum(jmt), :) , 'LineWidth', 3)
            plot( coords{jmt}(1,1), coords{jmt}(2,1), 'b*', 'MarkerSize', 10)
            LH( cNum(jmt) ) = plot(nan, nan, '-', 'Color', col( cNum(jmt), :), 'MarkerSize', 14, 'LineWidth', 10); L{ cNum(jmt) } = ['Line ', num2str( lineNum(jmt) )];    
        end
        hold off
        legend( LH, L)
        set( gcf, 'WindowState', 'maximized')
    
        if ~isempty(lineNum), nextIdx = lineNum(end)+1; else nextIdx=startIdx; end
    end
    % }}}

    % displayAllZplot {{{
    function displayAllZplot( im2, mts, plotflag)
    
    if ~plotflag, return, end

    numMT = length( mts);
    col = distinguishable_colors( numMT, {'w', 'k'});
    
    dispImg( im2.*imMask); set(gcf, 'NumberTitle', 'off', 'Name', 'allZ'); hold on;
    LH = []; L = {};;
    for jmt = 1 : numMT
        plot( mts(jmt).coords(1,:), mts(jmt).coords(2,:), 'LineWidth', 2, 'Color', col(jmt, :) ); hold on;
        LH( jmt) = plot(nan, nan, '-', 'Color', col( jmt, :), 'MarkerSize', 14, 'LineWidth', 10); L{ jmt} = ['Line', num2str( jmt)];    
    end
    title('allZ');
    legend( LH, L); hold off;

    end
    % }}}

    % analyzeLine {{{
    function mts = analyzeLine( mts, imMask )

        % fields: coords, coefX, coefY, CoM, orient, endPt, length, image
        % each line will be approximated with a 1st order parametric polynomial( linear line) in 2D
        nummt = length(mts);
        for j1 = 1 : nummt 
            x = mts(j1).coords(1,:);
            y = mts(j1).coords(2,:);
            t = linspace( 0, 1, length(x) );
            cx = polyfit( t, x, 1); cy = polyfit(t, y, 1);
            mts(j1).coefX = cx; mts(j1).coefY = cy;
            mts(j1).com = [ polyval(cy, 0.5), polyval(cx, 0.5) ];
            mts(j1).orient = atan( diff(polyval(cy, [0,1]))/diff(polyval(cx, [0,1])) );
            mts(j1).endPt = [polyval( cx, [0,1]) ; polyval(cy, [0,1]) ];
            mts(j1).length = norm( [cx(1), cy(1) ]);

            [cX, cY] = interpolateCoords( x, y, 10);
            indMT = sub2ind( size(imMask), round(cY), round(cX) );
            imLine = logical( 0*imMask);
            imLine( indMT) = 1;
            mts(j1).image = imLine;

        end
    
    end
    % }}}

    % findLinkCosts {{{
    function [costs, mts, imLinkArray] = findLinkCosts( mts, costAccept, plotflag)

        % Initialization {{{ 
        displayLinkInfo = 0;

        % Link Restrictions
        max_frame_dist = 1;
        nolink_removal = 1;

        % Link Parameters
        cost_accept = costAccept;
        axis_dist_max = 3;
        axis_phi_max = deg2rad( 22.5);
        coi_dist_max =  30;
        line_length_max = 100;
        line_length_min = 10;

        nummt = length(mts);
        costs = Inf( nummt);
        imLinkArray = cell( nummt, nummt);

        % }}}

        % Cost Calculation Functions: {{{
        function cost = cost_axis_dist(axis_dist, axis_dist_max, cost_accept, order)
            cost = ( cost_accept / axis_dist_max^order ) * axis_dist^order;
        end

        function cost = cost_axis_phi(axis_phi, axis_phi_max, cost_accept, order)
            cost = ( cost_accept / axis_phi_max^order) * axis_phi^order;
        end

        function cost = cost_coi_dist(coi_dist, coi_dist_max, cost_accept, order)
            cost = ( cost_accept / coi_dist_max^order) * coi_dist^order;
        end

        function cost = cost_line_length(line_length, line_length_max, line_length_min, cost_accept, order)
            cost = cost_accept - ( cost_accept / (line_length_max-line_length_min)^order) * line_length^order;
            if line_length < line_length_min, cost = cost_accept; end
            if line_length > line_length_max, cost = 0; end
        end
        % }}}

        % Loop over Filaments 
        for j1 = 1 : nummt, for j2 = 1: nummt
        
            % Conditional Linking {{{
            blockLink = 0;
            
            % Disable same-frame Linking
            if mts(j1).z == mts(j2).z
                blockLink = 1;
            end
            % Disable far-away Linking
            if abs( mts(j1).z - mts(j2).z ) > max_frame_dist
                blockLink = 1;
            end
            % Disable Linking with itself
            if j1 == j2
                blockLink = 1;
            end
            % Disable double Linking (j1 to j2, and j2 to j1) : eliminiates redundant information
            if j1 > j2
                blockLink = 1;
            end

            if blockLink, continue, end
            % }}} 
            
            % Same Frame {{{
            if mts(j1).z == mts(j2).z

                % find closest ends
                [pt, ptIdx, linkDist] = findClosestEnds( mts(j1), mts(j2) )
                
                % find orientation of lines and the link
                [phi_1, phi_2, phi_link] = findOrientationsOfLinesAndLink( ptIdx, pt, mts(j1), mts(j2))
                
                % Phi Cost
                phi_diff = abs( phi_link - phi_1) + abs(phi_link - phi_2);
                cost_phi = scale_phi(phi_diff)*( 5 - 4*linkDist/(2*max_link_length) );
                
                % Length Cost
                cost_dist = scale_dist_sameZ(linkDist);
                
                % Plotting {{{
                if plotflag
                     figure;
                     imagesc( max(imGauss, [], 3).*imMask); colormap gray; hold on;
                     plot( polyval( mts(j1).coefX, 0:0.05:1), polyval( mts(j1).coefY, 0:0.05:1), 'b-', 'linewidth', 3)
                     plot( polyval( mts(j2).coefX, 0:0.05:1), polyval( mts(j2).coefY, 0:0.05:1), 'r-', 'linewidth', 3)
                     plot( pt(1,:), pt(2,:), 'g-', 'linewidth', 3)
                     text( 20, 20, sprintf('phi1=%.2f \n phi2=%.2f \n phiLink=%.2f', phi_1, phi_2, phi_link), 'HorizontalAlignment', 'left', 'Color', 'White', 'FontSize', 14 ) 
                     text( size(imMask, 1), 20, sprintf('cost=%.2f \n cost axis=%.2f \n cost orient=%.2f', cost_dist+cost_orient, cost_dist, cost_orient), 'HorizontalAlignment', 'right', 'Color', 'White', 'FontSize', 14 ) 
                     set(gcf, 'NumberTitle', 'off', 'Name', sprintf('cost_match_%d_with_%d' , j2, j1))
                     hold off
                end
                % }}}

            end
            % }}}

            % Different Frame {{{
            if mts(j1).z ~= mts(j2).z
                
                % Measurements for Cost Calculation
                % 1. Axis Dist
                % 2. Axis Phi
                % 3. COI Dist
                % 4. Line Length

                % 1. Axis Dist
                % find distance of mt1 from axis of mt2 and vice versa. axis dist is the mean of the two.
                axis_dist = mean( [ findMinimumDistanceBetweenLines( mts(j1), mts(j2) ), findMinimumDistanceBetweenLines( mts(j2), mts(j1) ) ] );
                costAxisDist = cost_axis_dist( axis_dist, axis_dist_max, cost_accept, 1);

                % 2. Axis Phi 
                % find difference in orientations of axis of the 2 microtubules
                axis_phi = abs( mts(j1).orient - mts(j2).orient);
                costAxisPhi = cost_axis_phi( axis_phi, axis_phi_max, cost_accept, 1);
               
                % 3. COI Dist
                % find distance from COI_1 to COI_2
                coi_dist = norm( mts(j1).com - mts(j2).com );
                costCoiDist = cost_coi_dist( coi_dist, coi_dist_max, 0.5*cost_accept, 1);
                
                % 4. Line Length
                % find mean line length of the linking mts
                line_length = mean( [mts(j1).length, mts(j2).length] );
                costLineLength = cost_line_length( line_length, line_length_max, line_length_min, cost_accept, 2);
               
                if displayLinkInfo
                    disp( sprintf( 'link: %d to %d', j1, j2) )
                    disp( sprintf( 'axis_dist = %.2f', axis_dist) )
                    disp( sprintf( 'axis_phi= %.2f', axis_phi) )
                    disp( sprintf( 'coi_dist = %.2f', coi_dist) )
                    disp( sprintf( 'line_length = %.2f', line_length) )
                    disp( '----------------------------')
                end

                costs(j1, j2) = costAxisDist + (abs(costAccept-costCoiDist)/costAccept) *costAxisPhi + costCoiDist;
                
                % Link Image Creation {{{
                % find image of connection
                switch ( coi_dist - max( [mts(j1).length/2, mts(j2).length/2 ] ) ) < 0
                    case 1
                    % overlapping
                    imLink = logical( mts(j1).image + mts(j2).image);
                    case 0
                    % non-overlapping
                    [end1, end2, ~] = findClosestEnds( mts(j1), mts(j2)); end1 = round(end1); end2 = round(end2);
                    imLink = 0*mts(j1).image; imLink( end1(2), end1(1) ) = 1; imLink( end2(2), end2(1) ) = 1; imLink = bwconvhull( imLink);
                    imLink = logical(imLink + mts(j1).image + mts(j2).image);
                    otherwise
                    error('how is this possible?')
                end
                imLink = imdilate( imLink, strel('disk', 2) );
                imLinkArray{ j1, j2} = imLink; imLinkArray{ j2, j1} = imLink;
                %  }}}

            end
            % }}}

            % Plotting {{{
            if plotflag
                 figure;
                 subaxis(1, 2, 1)
                 imagesc( max(imGauss, [], 3).*imMask); colormap gray; hold on;

                 plot( polyval( mts(j1).coefX, 0:0.05:1), polyval( mts(j1).coefY, 0:0.05:1), 'b-', 'linewidth', 3)
                 plot( polyval( mts(j1).coefX, -2:0.05:2), polyval( mts(j1).coefY, -2:0.05:2), 'b.', 'linewidth', 3)
                 plot( polyval( mts(j2).coefX, 0:0.05:1), polyval( mts(j2).coefY, 0:0.05:1), 'r-', 'linewidth', 3)
                 plot( polyval( mts(j2).coefX, -2:0.05:2), polyval( mts(j2).coefY, -2:0.05:2), 'r.', 'linewidth', 3)

                 text( size(imMask, 1), 20, sprintf('cost = %.2f \n cost axis = %.2f \n cost phi = %.2f \n cost coi = %.2f \n cost length - %.2f', costs(j1, j2), costAxisDist, costAxisPhi, costCoiDist, costLineLength), 'HorizontalAlignment', 'right', 'Color', 'White', 'FontSize', 14 ) 
                 hold off

                 subaxis(1, 2, 2)
                 imagesc( imLink), colormap gray, hold on;
                 plot( mts(j1).coords(1,:), mts(j1).coords(2,:), 'b', 'linewidth', 3)
                 plot( mts(j2).coords(1,:), mts(j2).coords(2,:), 'r', 'linewidth', 3)
                 hold off;
                 set(gcf, 'NumberTitle', 'off', 'Name', sprintf('cost_match_%d_with_%d' , j1, j2), 'WindowState','maximized')
            end
            % }}}

        end, end
    
        % findClosestEnds {{{
        function [end1, end2, minDist] = findClosestEnds( mt1, mt2 )
            % there are 2 endpoints per line. we find the 2 closest endpoints
            minDist = Inf;
            end1 = mt1.endPt(:,1); end2 = mt2.endPt(:,1);
            for p1 = 1:2, for p2 = 1:2
                dist = norm( mt1.endPt(:, p1) - mt2.endPt(:, p2) );
                if dist < minDist
                    minDist = dist;
                    end1 = mt1.endPt(:,p1); end2 = mt2.endPt(:, p2);
                end
            end, end
        end
        % }}}
                
        % findOrientationsOfLinesAndLink {{{
        function [phi_1, phi_2, phi_link] = findOrientationsOfLinesAndLink( ptIdx, pt, mt1, mt2)
            % orientation of the mts and the link
            % phi of 1
            if ptIdx(1) == 1
                riseRun = mt1.endPt( :, 1 ) - mt1.endPt(:,2); 
            elseif ptIdx(1) == 2
                riseRun = mt1.endPt( :, 2 ) - mt1.endPt(:,1); 
            end
            phi_1 = atan2( riseRun(2),riseRun(1) );

            % phi of 2
            if ptIdx(2) == 1
                riseRun = mt2.endPt( :, 2 ) - mt2.endPt(:,1); 
            elseif ptIdx(2) == 2
                riseRun = mt2.endPt( :, 1 ) - mt2.endPt(:,2); 
            end
            phi_2 = atan2( riseRun(2),riseRun(1) );

            % find Orientation of the Link
            phi_link = atan2( diff( pt(2,:) ), diff(pt(1,:) ) );

        end
        % }}}

        % findMinimumDistanceBetweenLines {{{
        function axisDist = findMinimumDistanceBetweenLines( mt1, mt2)

            % pick a point on mt_j2 and find its minimum perp distance from the axis of mt_j1
            numTestPts = 15;
            s2 = linspace( 0,1, numTestPts);
            x2 = polyval( mt2.coefX, s2); y2 = polyval( mt2.coefY, s2);
            dists = [];
            for jp = 1 : numTestPts 
                % find corresponding point on the axis of mt_j1 which is closest to this point
                s1 = ( mt1.coefX(1) * ( x2(jp) - mt1.coefX(2) ) + mt1.coefY(1) * ( y2(jp) - mt1.coefY(2) ) ) / ( mt1.coefX(1)^2 + mt1.coefY(1)^2 );
                x1 = polyval( mt1.coefX, s1); y1 = polyval( mt1.coefY, s1);
                dists = [dists, norm( [x1,y1] - [x2(jp), y2(jp)]) ];
            end
            
            % axis distance
            axisDist = min(dists);

        end
        % }}}

    end
    % }}}

    % linkMts {{{
    function [filamentBank, filamentID, mts] = linkMts( mts, costMat, costAccept, imMask, linkArray, plotflag)
        
        filamentBank = {};
        filamentID = {};
        nolinkRemoval = 1;
        nummt = length(costMat);
        
        % Assemble New Filaments {{{
        % Assign ID to each filament.
        for jmt = 1:nummt
            mts(jmt).id = jmt;
            mts(jmt).imageNew = mts(jmt).image;
        end
        % if a filament is to be linked, we'll assign the same id to the linking couple.
        for jR = 1:nummt
            idxLink = find(costMat(jR,:) < costAccept);
            % find link label
            linklabel = min( [mts([jR, idxLink]).id]);
            % find link Image
            linkImage = mts(jR).imageNew;
            for jLink = idxLink
                linkImage = logical( linkImage + linkArray{ jR, jLink} );
            end
            % assign label and image
            for jLink = [jR,idxLink]
                mts(jLink).id = linklabel;
            end
            matches = find( [mts.id] == linklabel);
            for jLink = matches
                mts(jLink).imageNew = linkImage;
            end
%             figure; imagesc( cat(3, mts(jR).imageNew, linkImage, 0*linkImage) ), title( sprintf('mt %d to mt %d', jR, idxLink) )
        end
        % }}}

        disp( [mts.id])

        % Finalize New Filaments {{{
        % Find the unique objects
        uniqueID = unique( [mts.id]);
        for jmt = 1 : length(uniqueID)
            cid = uniqueID(jmt);
            idxMT = find( [mts.id]==cid);
            if length( idxMT) == 1 && nolinkRemoval
                mts(idxMT).remove = 1;
            else
                filamentID = { filamentID{:}, idxMT};
                filamentBank = {filamentBank{:}, mts( idxMT(1) ).imageNew}; 
            end
        end
        numFilaments = length( filamentBank);
        % }}}
        
        % Plotting {{{
        if plotflag
            dispImg( filamentBank{:}, [1 length(filamentBank)]);
            set( gcf, 'NumberTitle', 'off', 'Name', 'Linked Filaments: Images');

            figure;
            for jmt = 1: numFilaments
                subaxis(1, numFilament, jmt);
                imagesc( max(imGauss, [], 3) ); colormap gray; axis equal;
                nCol = distinguishable_colors( length( filamentID{jmt} ), {'w', 'k'} );
                for jid = filamentID{ jmt};
                    plot( mts(jid).coords(1,:), mts(jid).coords(2,:), 'Color', nCol(:,jid), 'LineWidth', 2); 
                end
            end
            set( gcf, 'NumberTitle', 'off', 'Name', 'Linked Filaments');
        end
        % }}} 

    end
    % }}}
   
    % findFinalMT {{{
    function [coords, orients] = findFinalMT( imagesMT, imagesID, mts, imGauss, minMTlength, plotflag)

        % Initialization {{{
        coords = {};
        orients = [];
        imMIP = max(imGauss, [], 3);
        bkg_thresholding = 0;
        nummt = length(imagesMT);
        plotflags.success=0;
        plotflags.fail=0;
        % }}}
       
       % Plotting {{{
        if plotflag
            dispImg( imagesMT{:}, [1 nummt])
            for jmt = 1: nummt
                subplot(1, nummt, jmt)
                hold on;
                numid = length( imagesID{jmt} );
                col = distinguishable_colors( numid, {'w', 'k'});
                for jid = 1: numid
                    cid = imagesID{jmt}(jid);
                    plot( mts( cid).coords(1,:), mts(cid).coords(2,:), 'Color', col(jid, :), 'LineWidth', 3);
                end
            end
            set(gcf, 'NumberTitle', 'off', 'Name', 'Final Linked Filaments')
        end
        % }}} 

        % Now we can feed each of these images and detect microtubule along these profiles
        for jim = 1 : nummt 
            
            % Find MT {{{
            try
            imMT = imagesMT{ jim}.*imMIP; 
            catch, size(imagesMT{jim}), size(imMIP), end 

            %  find mtoc for this microtubule
            [ intMTOC, idxMTOC] = max( imMT(:) ); [ymtoc, xmtoc] = ind2sub( size(imMT), idxMTOC);

            [~, success, orientStart] = estimateNextPoint( [ xmtoc; ymtoc], 0, 5, 15, 540, imMT, bkg_thresholding) ;

            if ~success
                error('I"m broken...oh wont you fix me?')
            end

            phi= [orientStart, orientStart+pi];
            for jAng = 1:2
                cMT = estimateMicrotubulePoints( [xmtoc; ymtoc], phi(jAng), imMT, 5, 20, 70, bkg_thresholding, plotflags);
                if size( cMT, 2)*step > minMTlength
                    coords = { coords{:}, cMT};
                    orients = [ orients, phi(jAng)];
                end
            end
            % }}}

        end

    end
    % }}}


end
