classdef Methods

    properties
    end

    methods ( Static = true)
        
        % FindCurves {{{ 
        function coords = FindCurves( imageIn, varargin)

            disp('Searching for Curves')

            % Parse Arguments
            opts = parseArgs( imageIn, varargin{:});

            % Dimension
            dim = length( size( imageIn));

            % Apply mask to image
        %     imMasked = imageIn .* opts.Mask;

            % Get Filtered and Steerable images
            [imGauss, imSteer] = Methods.FilterImageForCurveDetection( imageIn, 'Mask', opts.Mask);

            % Find dimension specific curves
            vars = Methods.Opts2VarsExcept( opts, 'ImageSteer');
            if dim == 3
                coords = Methods.FindCurves3D( imGauss, 'ImageSteer', imSteer, vars{:});
            elseif dim == 2
                coords = Methods.FindCurves2D( imGauss, 'ImageSteer', imSteer, vars{:});
            end

            % parseArgs {{{
            function opts = parseArgs( image, varargin)
                % Possible Arguments
                % Mask
                % Visibility
                % FieldOfView
                % StepSize
                % MinLength
                % MaxCurves
                % CostAcceptLink
                % MaxFails
                % Plot
                % Verbose

                % Dimension of image
                dim = length( size( image) );

                % Default Values
                defaultMask = image ~= 0;
                defaultVisibility = 10;
                defaultFieldOfView = 40;
                defaultStepSize = 5;
                defaultMinLength = 10;
                defaultMaxCurves = 5;
                defaultCostAcceptLink = 10;
                defaultMaxFails = 2;
                defaultPlot = 0; % 0: no plots, 1: important plots, 2: debug-all plots
                defaultVerbose = 0;

                % Input Parser
                p = inputParser;
                p.KeepUnmatched=true;

                % Mask 
                validMask = @(x) length(size(x)) == dim && all( size(x) == size(image)) ;
                addParameter( p, 'Mask', defaultMask, validMask);

                % Function for positive integer
                validPosInt = @(x) length(x)==1 && isnumeric(x) && x>0; 

                % Visibility 
                addParameter( p, 'Visibility', defaultVisibility, validPosInt);
                % Field of View 
                addParameter( p, 'FieldOfView', defaultFieldOfView, validPosInt);
                % StepSize
                addParameter( p, 'StepSize', defaultStepSize, validPosInt);
                % MinLength
                addParameter( p, 'MinLength', defaultMinLength, validPosInt);
                % MaxCurves
                addParameter( p, 'MaxCurves', defaultMaxCurves, validPosInt);
                % CostAcceptLink
                addParameter( p, 'CostAcceptLink', defaultCostAcceptLink, validPosInt);
                % MaxFails
                addParameter( p, 'MaxFails', defaultMaxFails, validPosInt);
                % Plot 
                validPlot = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1 || x==2 || x==3); 
                addParameter( p, 'Plot', defaultPlot, validPlot);
                % Verbose 
                validVerbose = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1); 
                addParameter( p, 'Verbose', defaultVerbose, validVerbose);

                parse( p, varargin{:} );
                opts = p.Results;

            end
            % }}}

        end
        % }}}

        % FindCurves3D {{{
        function coords = FindCurves3D( imageIn, varargin)

            % To estimate the lines, we'll take the following approach:
            % 1) Start at the brightest point in the image, and find the optimal orientation of any underlying lines by radially integrating intensity.
            % 2) Once we have a start point and a direction, we continue searching for the underlying line by looking for intensity peaks until a threshold is met.
            % 3) We take this estimated microtubule, and the pixel region it spans and set it to the median value of the original image. We say that we have removed the intensity due to this microtubule.
            % 4) We repeat steps 1-3 until the microtubule estimates are smaller than a certain length
            % 5) We will fit a polynomial through each initial estimate (quadratic or
            %    cubic). This will be our final estimate for the microtubules.

            % Parse Arguments
            opts = parseArgs( imageIn, varargin{:});

            % Get Steerable image if not present
            if all( isnan( opts.ImageSteer) )
                [imageIn, opts.ImageSteer] = Methods.FilterImageForCurveDetection( imageIn, 'Mask', opts.Mask);
            end

            % Get image stats
            stats = Methods.GetImageStats( imageIn, opts.Verbose);

            cLine = 1; 
            for jZ = 1 : size( imageIn, 3)

                %  Z-specific params
                disp( sprintf('  Z-Slice = %d', jZ) )
                imGaussZ = imageIn(:,:,jZ);
                imSteerZ = opts.ImageSteer(:,:,jZ);

                % Find 2D curves in all Z frames
                vars = Methods.Opts2VarsExcept( opts, 'ImageSteer', 'ThreshInt', 'Mask', 'Verbose');
                [coords, out] = Methods.FindCurves2D( imGaussZ, 'ImageSteer', imSteerZ, 'ThreshInt', stats.ThreshLow, 'Verbose', 1, vars{:} );

                displayZplot( imGaussZ, imSteerZ, coords, cLine, opts.Plot>=2);
                if opts.Plot>=2, set( gcf, 'NumberTitle', 'off', 'Name', sprintf('z=%d', jZ) ), end

                % Store information about found curves in a useful format
                for jC = 1 : length(coords)
                    curves( cLine).coords = coords{jC};
                    curves( cLine).z = jZ;
                    curves( cLine).int = out.SeedInt(jC);
                    cLine = cLine + 1;
                end

            end

            displayAllZplot( max(imageIn, [], 3), curves, opts.Plot>=2);
                        
            % Analyze each detected line and save its information
            curves = analyzeLine( curves, opts.Mask);

            % find link costs between microtubules in structure
            [costs, curves, linkArray] = findLinkCosts( curves, opts.CostAcceptLink, opts.Plot>=3, imageIn);
            
            % link mts if the cost is below a certain threshold
            [imagesMT, imagesID, curves] = linkMts( curves, costs, opts.CostAcceptLink, imageIn, linkArray, opts.Plot>=3);
            
            
            % Find final curve
            [coords, ~] = findFinalMT( imagesMT, imagesID, curves, imageIn, opts.MinLength, opts.StepSize, opts.Plot>=1);
            
            % Add z coordinate
            for jseed = 1 : length( coords)
                
                % get seed (x,y) coordinates
                coord2 = coords{jseed}{1}(:,1);
                % get z coordinate at seed and duplicate it for all curve
                % coordinates
                [~, zidx] = max( imageIn, [],3);
                for jc = 1 : length(  coords{jseed})
                    coords{jseed}{jc}(3,:) = zidx(coord2(2), coord2(1));
                end
            end

            % data for paper figure
            curveFigData = 0;
            if curveFigData == 1
                sP = [pwd, filesep, 'PaperFigures/Estimation/'];
                save( [sP, 'mts&imgs.mat'], 'curves', 'imagesMT')
            end

            % parseArgs {{{
            function opts = parseArgs( image, varargin)
                % Possible Arguments
                % Mask
                % Visibility
                % FieldOfView
                % StepSize
                % MinLength
                % MaxCurves
                % CostAcceptLink
                % MaxFails
                % Plot
                % Verbose
                % ImageSteer
                
                % Dimension of image
                dim = length( size( image) );

                % Default Values
                defaultMask = image ~= 0;
                defaultVisibility = 10;
                defaultFieldOfView = 40;
                defaultStepSize = 5;
                defaultMinLength = 10;
                defaultMaxCurves = 5;
                defaultCostAcceptLink = 10;
                defaultMaxFails = 2;
                defaultPlot = 0; % 0: no plots, 1: important plots, 2: debug-all plots
                defaultVerbose = 0;
                defaultImageSteer = nan( size(image) );

                % Input Parser
                p = inputParser;
                p.KeepUnmatched=true;

                % Mask 
                validMask = @(x) length(size(x)) == dim && all( size(x) == size(image));
                addParameter( p, 'Mask', defaultMask, validMask);

                % Function for positive integer
                validPosInt = @(x) length(x)==1 && isnumeric(x) && x>0; 

                % Visibility 
                addParameter( p, 'Visibility', defaultVisibility, validPosInt);
                % Field of View 
                addParameter( p, 'FieldOfView', defaultFieldOfView, validPosInt);
                % StepSize
                addParameter( p, 'StepSize', defaultStepSize, validPosInt);
                % MinLength
                addParameter( p, 'MinLength', defaultMinLength, validPosInt);
                % MaxCurves
                addParameter( p, 'MaxCurves', defaultMaxCurves, validPosInt);
                % CostAcceptLink
                addParameter( p, 'CostAcceptLink', defaultCostAcceptLink, validPosInt);
                % MaxFails
                addParameter( p, 'MaxFails', defaultMaxFails, validPosInt);
                % Plot 
                validPlot = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1 || x==2|| x==3); 
                addParameter( p, 'Plot', defaultPlot, validPlot);
                % Verbose 
                validVerbose = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1); 
                addParameter( p, 'Verbose', defaultVerbose, validVerbose);

                % ImageSteer 
                validImageSteer = @(x) length(size(x)) == dim && all( size(x) == size(image));
                addParameter( p, 'ImageSteer', defaultImageSteer, validImageSteer);

                parse( p, varargin{:} );
                opts = p.Results;

            end
            % }}}

            % displayZplot {{{
            function displayZplot( im2D, imLines, coords, startIdx, plotflag)
                
                if ~plotflag
                    return
                end

                % Construct line numbers
                lineNum = startIdx : startIdx+length( coords)-1;
                lineNumIdx = 1: length(lineNum);
                if isempty( lineNum)
                    skipCurvePlot = 1;
                else, skipCurvePlot = 0; end
                
                LH=[]; L={};
                fS = size(im2D, 1);
                figure;
                subaxis(1,3,1); imagesc( im2D ); colormap gray; axis equal; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
                subaxis(1,3,2); imagesc( imLines ); colormap gray; axis equal; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
                subaxis(1,3,3); imagesc( im2D ); colormap gray; axis equal; hold on; set( gca, 'xlim', [1 fS], 'ylim', [1 fS] );
                
                if ~skipCurvePlot
                    col = distinguishable_colors( lineNumIdx(end), {'k', 'w'});
                end
                for jmt = 1 : length(coords)
                    plot( coords{jmt}(1,:), coords{jmt}(2,:), 'color', col( lineNumIdx(jmt), :) , 'LineWidth', 3)
                    plot( coords{jmt}(1,1), coords{jmt}(2,1), 'b*', 'MarkerSize', 10)
                    LH( lineNumIdx(jmt) ) = plot(nan, nan, '-', 'Color', col( lineNumIdx(jmt), :), 'MarkerSize', 14, 'LineWidth', 10); L{ lineNumIdx(jmt) } = ['Line ', num2str( lineNum(jmt) )];    
                end
                hold off
                if ~skipCurvePlot
                    legend( LH, L)
                end
                set( gcf, 'WindowState', 'maximized')
            
            end
            % }}}
                    
            % displayAllZplot {{{
            function displayAllZplot( im2, mts, plotflag)
            
                if ~plotflag, 
                    return, 
                end

                numMT = length( mts);
                col = distinguishable_colors( numMT, {'w', 'k'});
                
                dispImg( im2); set(gcf, 'NumberTitle', 'off', 'Name', 'allZ'); hold on;
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
                imMask = imMask(:,:,1);
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

                    [cX, cY] = Methods.InterpolateCoords( x, y, 10);
                    indMT = sub2ind( size(imMask), round(cY), round(cX) );
                    imLine = logical( 0*imMask);
                    imLine( indMT) = 1;
                    mts(j1).image = imLine;

                end

            end
            % }}}

            % Linking {{{
            % findLinkCosts {{{
            function [costs, mts, imLinkArray] = findLinkCosts( mts, costAccept, plotflag, imGauss)

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
%                     if mts(j1).z == mts(j2).z
%                         blockLink = 1;
%                     end
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
                        
                        % Measurements for Cost Calculation
                        % 1. Axis Dist
                        % 2. Axis Phi

                        % 1. Axis Dist
                        % find distance of mt1 from axis of mt2 and vice versa. axis dist is the mean of the two.
                        axis_dist = mean( [ findMinimumDistanceBetweenLines( mts(j1), mts(j2) ), findMinimumDistanceBetweenLines( mts(j2), mts(j1) ) ] );
                        costAxisDist = cost_axis_dist( axis_dist, axis_dist_max, cost_accept, 1);

                        % 2. Axis Phi 
                        % find difference in orientations of axis of the 2 microtubules
                        axis_phi = abs( mts(j1).orient - mts(j2).orient);
                        axis_phi = min( [pi - axis_phi, axis_phi, pi + axis_phi]);
                        costAxisPhi = cost_axis_phi( axis_phi, axis_phi_max, cost_accept, 1);

                        if displayLinkInfo
                            disp( sprintf( 'link: %d to %d', j1, j2) )
                            disp( sprintf( 'axis_dist = %.2f', axis_dist) )
                            disp( sprintf( 'axis_phi= %.2f', axis_phi) )
                            disp( '----------------------------')
                        end

                        costs(j1, j2) = costAxisDist + costAxisPhi;

                        % find image of connection
                        [end1, end2, ~] = findClosestEnds( mts(j1), mts(j2)); end1 = round(end1); end2 = round(end2);
                        imLink = 0*mts(j1).image; imLink( end1(2), end1(1) ) = 1; imLink( end2(2), end2(1) ) = 1; imLink = bwconvhull( imLink);
                        imLink = logical(imLink + mts(j1).image + mts(j2).image);
                        imLink = imdilate( imLink, strel('disk', 2) );
                        imLinkArray{ j1, j2} = imLink; imLinkArray{ j2, j1} = imLink;

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
                        axis_phi = min( [pi - axis_phi, axis_phi, pi + axis_phi]);
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
                         imagesc( max(imGauss, [], 3)); colormap gray; hold on;

                         plot( polyval( mts(j1).coefX, 0:0.05:1), polyval( mts(j1).coefY, 0:0.05:1), 'b-', 'linewidth', 3)
                         plot( polyval( mts(j1).coefX, -2:0.05:2), polyval( mts(j1).coefY, -2:0.05:2), 'b.', 'linewidth', 3)
                         plot( polyval( mts(j2).coefX, 0:0.05:1), polyval( mts(j2).coefY, 0:0.05:1), 'r-', 'linewidth', 3)
                         plot( polyval( mts(j2).coefX, -2:0.05:2), polyval( mts(j2).coefY, -2:0.05:2), 'r.', 'linewidth', 3)

%                          text( size(imGauss, 1), 20, sprintf('cost = %.2f \n cost axis = %.2f \n cost phi = %.2f \n cost coi = %.2f \n cost length - %.2f', costs(j1, j2), costAxisDist, costAxisPhi, costCoiDist, costLineLength), 'HorizontalAlignment', 'right', 'Color', 'White', 'FontSize', 14 ) 
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
            function [filamentBank, filamentID, mts] = linkMts( mts, costMat, costAccept, imGauss, linkArray, plotflag)

                filamentBank = {};
                filamentID = {};
                keepUnlinkedIfGood = 1;
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

%                 disp( [mts.id])
                stats = Methods.GetImageStats( imGauss);

                % Finalize New Filaments {{{
                % Find the unique objects
                uniqueID = unique( [mts.id]);
                for jmt = 1 : length(uniqueID)
                    cid = uniqueID(jmt);
                    idxMT = find( [mts.id]==cid);
                    
                    % If unlinked
                    if length( idxMT) == 1 
                        
                        % Keep if curve is long and bright and on an edge
                        % frame
                        if keepUnlinkedIfGood ...
                                && mts(idxMT).length > 30 ...
                                && mts(idxMT).int > stats.ThreshHigh ...
                                && ( mts(idxMT).z == 1 || mts(idxMT).z == size(imGauss,3) )
                            filamentID = { filamentID{:}, idxMT};
                            filamentBank = {filamentBank{:}, mts( idxMT(1) ).imageNew}; 
                        else % remove
                            mts(idxMT).remove = 1;
                        end
                        
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
                        subaxis(1, numFilaments, jmt);
                        imagesc( max(imGauss, [], 3) ); colormap gray; axis equal; hold on;
                        nCol = distinguishable_colors( length( filamentID{jmt} ), {'w', 'k'} );
                        for jid = 1: length( filamentID{ jmt} )
                            id = filamentID{jmt}(jid);
                            plot( mts(id).coords(1,:), mts(id).coords(2,:), 'Color', nCol(jid,:), 'LineWidth', 2); 
                        end
                        hold off;
                    end
                    set( gcf, 'NumberTitle', 'off', 'Name', 'Linked Filaments');
                end
                % }}} 

            end
            % }}}

            % findFinalMT {{{
            function [coords, orients] = findFinalMT( imagesMT, imagesID, mts, imGauss, minMTlength, step, plotflag)

                % Initialization {{{
                coords = {};
                orients = {};
                imMIP = max(imGauss, [], 3);
                bkg_thresholding = 0;
                nummt = length(imagesMT);
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
                    coordsAster = {};
                    orientsAster = [];
                    for jAng = 1:2
                        cMT = estimateMicrotubulePoints( [xmtoc; ymtoc], phi(jAng), imMT, 5, 20, 70, bkg_thresholding);
                        if size( cMT, 2)*step > minMTlength
                            coordsAster = { coordsAster{:}, cMT};
                            orientsAster = [ orientsAster, phi(jAng)];
                            % coords = { coords{:}, cMT};
                            % orients = [ orients, phi(jAng)];
                        end
                    end
                    coords = { coords{:}, coordsAster};
                    orients = { orients{:}, orientsAster};
                    % }}}

                end
                
                if plotflag
                    dispImg( imMIP); hold on;
                    col = distinguishable_colors( length(coords), {'w', 'k', 'r'});
                    for jb = 1 : length(coords)
                        for jmt = 1: length( coords{jb})
                            plot( coords{jb}{jmt}(1,:), coords{jb}{jmt}(2,:), 'Color', col(jb, :), 'LineWidth', 3);
                            plot( coords{jb}{jmt}(1,1), coords{jb}{jmt}(2,1), 'r*', 'LineWidth', 5, 'MarkerSize', 10);
                        end
                    end
                    set(gcf, 'NumberTitle', 'off', 'Name', 'Final Curves')
                end

            end
            % }}}
            % }}}

        end
        % }}}

        % FindCurves2D {{{
        function [coords, outout] = FindCurves2D( imageIn, varargin)

            % Parse Arguments
            opts = parseArgs( imageIn, varargin{:});

            % Get Steerable image if not present
            if all( isnan( opts.ImageSteer) )
                [imageInFilt, opts.ImageSteer] = Methods.FilterImageForCurveDetection( imageIn, 'Mask', opts.Mask);
            end

            % Gaussian Steerable Image
%             imSteerGauss = imSteer.*imageIn;

            % Get image stats
            stats = Methods.GetImageStats( imageIn, opts.Verbose);

            % Initialize counters and arrays
            coords = {};
            fails = 0;
            iter = 1;
            intSeed = 1.2*opts.ThreshInt;
            imSearch = opts.ImageSteer;
            outout.Seed = [];
            outout.SeedInt = [];

            % Iteratively search for curves
            while (intSeed > opts.ThreshInt) && iter <= opts.MaxCurves && fails <= opts.MaxFails 

                vars = Methods.Opts2VarsExcept( opts);
                [coord, ~, out] = Methods.FindCurve2D( imSearch, vars{:});
                intSeed = out.IntSeed;

                % If Curve found
                if out.success

                    % Find the real seed 
                    [cX,cY] = Methods.InterpolateCoords( coord(1,:), coord(2,:), ceil(out.Length/size(coord,2)) ); 
                    indXY = sub2ind( size(imageIn), round(cY), round(cX) );
                    
                    % Trim curve from both ends to ensure curve intensity requirements are
                    % met
                    ind0 = find( imageIn(indXY) > stats.ThreshLow, 1, 'first');
                    ind1 = find( imageIn(indXY) > stats.ThreshLow, 1, 'last');
                    
                    if isempty(ind0) && isempty(ind1)
                        % Seed location and intensity
                        [intSeed, IdxSeed] = max( imageIn( indXY ));
                        [ Seed(2), Seed(1)] = ind2sub( size(imageIn),indXY(IdxSeed) );
                        out.Length = 0;
                    else
                        coord = [];
                        coord(1,:) = cX( ind0:ind1);
                        coord(2,:) = cY( ind0:ind1);
                        indXY = indXY(ind0:ind1);
                    
                        % Seed location and intensity
                        [intSeed, IdxSeed] = max( imageIn( indXY ));
                        [ Seed(2), Seed(1)] = ind2sub( size(imageIn),indXY(IdxSeed) );
                        out.Length = sum( sqrt( diff( coord(1,:)).^2 + diff( coord(2,:)).^2 ) );

                    end
                    
                    % Remove curve and seed from image 
                    imRm = Methods.RemoveSeedFromImage( imSearch, Seed, 3, opts.Mask);
                    imRm = Methods.RemoveCurveFromImage( imRm, coord, opts.Mask);

                    if opts.Verbose
                        disp( sprintf( '   iter = %d, int_Seed = %.2f, int_Curve = %.2f' , iter, intSeed, mean(imageIn(indXY)) ) )
                    end
                    
                    % display iteration plot if flag activated
                    displayIterationPlot( opts.ImageSteer, imSearch, imRm, coord, out.Seed, iter, opts.Plot>=2); 

                    % Accept curve if length bigger than some threshold and intensity is biggger than threshold
                    if out.Length>opts.MinLength && intSeed > opts.ThreshInt 
                        coords = {coords{:}, coord};
                        outout.Seed = [ outout.Seed ; Seed];
                        outout.SeedInt = [outout.SeedInt ; intSeed];
                    end
                    imSearch = imRm;

                else
                    
                    % Find bright seed
                    [ ~, idxMax] = max( imageIn(:) );
                    [ yMax, xMax] = ind2sub( size(imageIn), idxMax);
                    Seed = [xMax; yMax];
                    
                    % Remove seed from image
                    imSearch = Methods.RemoveSeedFromImage( imSearch, Seed, 3, opts.Mask);

                    % Update fails
                    fails = fails+1;
                end

                
                if opts.Verbose
                    disp( sprintf( '   Iter = %d , Fails = %d', iter, fails) )
                end
                iter = iter +1;

            end

            outout.coords = coords;           

            % parseArgs {{{
            function opts = parseArgs( image, varargin)
                % Possible Arguments
                % Mask
                % Visibility
                % FieldOfView
                % StepSize
                % MinLength
                % MaxCurves
                % MaxFails
                % Plot
                % Verbose
                % ImageSteer
                % ThreshInt
                
                % Dimension of image
                dim = length( size( image) );

                % Default Values
                defaultMask = image ~= 0;
                defaultVisibility = 10;
                defaultFieldOfView = 40;
                defaultStepSize = 5;
                defaultMinLength = 10;
                defaultMaxCurves = 5;
                defaultMaxFails = 2;
                defaultPlot = 0; % 0: no plots, 1: important plots, 2: debug-all plots
                defaultVerbose = 0;
                defaultImageSteer = zeros( size(image) );

                % default ThreshInt
                imVals = image( image ~= 0); imVals = imVals(:);
                defaultThreshInt = mean(imVals) + 3*std(imVals);

                % Input Parser
                p = inputParser;
                p.KeepUnmatched=true;

                % Mask 
                validMask = @(x) length(size(x)) == dim && all( size(x) == size(image));
                addParameter( p, 'Mask', defaultMask, validMask);

                % Function for positive integer
                validPosInt = @(x) length(x)==1 && isnumeric(x) && x>0; 

                % Visibility 
                addParameter( p, 'Visibility', defaultVisibility, validPosInt);
                % Field of View 
                addParameter( p, 'FieldOfView', defaultFieldOfView, validPosInt);
                % StepSize
                addParameter( p, 'StepSize', defaultStepSize, validPosInt);
                % MinLength
                addParameter( p, 'MinLength', defaultMinLength, validPosInt);
                % MaxCurves
                addParameter( p, 'MaxCurves', defaultMaxCurves, validPosInt);
                % MaxFails
                addParameter( p, 'MaxFails', defaultMaxFails, validPosInt);
                % Plot 
                validPlot = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1 || x==2|| x==3); 
                addParameter( p, 'Plot', defaultPlot, validPlot);
                % Verbose 
                validVerbose = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1) ;
                addParameter( p, 'Verbose', defaultVerbose, validVerbose);

                % ImageSteer 
                validImageSteer = @(x) ( length(size(x)) == dim) && ( all( size(x) == size(image)) );
                addParameter( p, 'ImageSteer', defaultImageSteer, validImageSteer);

                % ThreshInt 
                addParameter( p, 'ThreshInt', defaultThreshInt, validPosInt);

                parse( p, varargin{:} );
                opts = p.Results;

            end
            % }}}

            % displayIterationPlot {{{
            function displayIterationPlot( imOriginal, imPrior, imPost, coordMT, mtoc, iter, plotflag); 

                if ~plotflag,
                    return
                end

                figure('NumberTitle', 'off', 'Name', ['Iter_' num2str(iter)] )
                subplot(221)
                imagesc(imOriginal); colormap gray; title('Original Image')
                subplot(222)
                imagesc( imPrior); colormap gray; hold on; plot( mtoc(1), mtoc(2), 'b*', 'MarkerSize', 14); hold off;
                title( 'Search Image: Seed' )
                subplot(223)
                imagesc( imPrior); colormap gray; hold on; plot( coordMT(1, :), coordMT(2,:), 'r-', 'LineWidth', 2); plot( mtoc(1), mtoc(2), 'b*', 'MarkerSize', 14); hold off
                title( 'Curve found' )
                subplot(224)
                imagesc(imPost); colormap gray; title('Image for next iteration')

            end
            % }}}

        end
        % }}}

        % FindCurve2D {{{
        function [coordEndEnd, coordFromSeed, out] = FindCurve2D( imageIn, varargin)
            % FindCurve2D: finds a single curve given an intensity image
            % returns end end curve coords, and double coords from a possible seed

            % Parse Arguments
            opts = parseArgs( imageIn, varargin{:});

            % How to find a Curve:
            %   Find brightest spot in image
            %   Check for passing curves 
            %   If curve exists, then find it on both sides of bright spot
            %   Get coordinates for end to end curve

            % Find Bright Point in Image
            imG = imgaussfilt( imageIn, 1);
%             [ out.IntSeed, idxMax] = max( imageIn(:) );
            [ out.IntSeed, idxMax] = max( imG(:) );
            [ yMax, xMax] = ind2sub( size(imageIn), idxMax);
            out.Seed = [xMax; yMax];

            % Does curve  pass through this Bright Point? 
%             pf.success = 1; pf.fail = 1;
            [~, out.success, orientStart] = Methods.estimateNextCurveCoord( out.Seed, 0, opts.StepSize, opts.Visibility, 540, imG);
           
            out.Length = NaN;
            coordEndEnd = [];
            coordFromSeed = cell(1,2);
            % if curve passes, search for curve on both sides 
            if out.success
                
                % Set the two opposite angles for search
                angOrient = [ orientStart, orientStart+pi]; 
                for jAng = 1:2
                    coordFromSeed{jAng} = Methods.estimateCurveCoords( out.Seed, angOrient( jAng), imageIn, opts.StepSize, opts.Visibility, opts.FieldOfView, 1);
                end

                % Treat as single line and get ordered end to end coordinates
                coordEndEnd = [ flip( coordFromSeed{1}, 2), coordFromSeed{2} ]; 
                coordEndEnd( :, size(coordFromSeed{1}, 2) ) = [];

                % Length of end-end curve 
                out.Length = (size(coordEndEnd, 2)-1) * opts.StepSize;

            end

            % parseArgs {{{
            function opts = parseArgs( image, varargin)
                % Possible Arguments
                % Visibility
                % FieldOfView
                % StepSize
                % MinLength
                % Plot
                % Verbose
                
                % Dimension of image
                dim = length( size( image) );

                % Default Values
                defaultVisibility = 10;
                defaultFieldOfView = 40;
                defaultStepSize = 5;
                defaultMinLength = 10;
                defaultPlot = 0; % 0: no plots, 1: important plots, 2: debug-all plots
                defaultVerbose = 0;

                % Input Parser
                p = inputParser;
                p.KeepUnmatched=true;

                % Function for positive integer
                validPosInt = @(x) length(x)==1 && isnumeric(x) && x>0; 

                % Visibility 
                addParameter( p, 'Visibility', defaultVisibility, validPosInt);
                % Field of View 
                addParameter( p, 'FieldOfView', defaultFieldOfView, validPosInt);
                % StepSize
                addParameter( p, 'StepSize', defaultStepSize, validPosInt);
                % MinLength
                addParameter( p, 'MinLength', defaultMinLength, validPosInt);
                % Plot 
                validPlot = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1 || x==2|| x==3); 
                addParameter( p, 'Plot', defaultPlot, validPlot);
                % Verbose 
                validVerbose = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1) ;
                addParameter( p, 'Verbose', defaultVerbose, validVerbose);

                parse( p, varargin{:} );
                opts = p.Results;

            end
            % }}}

        end
        % }}}

        % estimateCurveCoords {{{
        function coordsMT  = estimateCurveCoords( startPt, startOrient, imageMT, stepSize, visibility, fieldOfVision, bkg_thresholding, plotflags)
            % EstimateMicrotubulePoints: estimates points along
            % microtubule.
            % We take the following approach:
            % 1) Nucleate at the MTOC. We find the angle that gives the maximum in a
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

            if nargin < 7, bkg_thresholding = 1; end
            if nargin < 8, plotflags.success = 0; plotflags.fail=0; end
            
            % clear out any previously stored estimated points
             
            success = 1; iter = 1;
            max_mt_length = 100; % max MT length approx
            max_iter = round(max_mt_length/stepSize);
            orientationOld = startOrient;
            coordCurr = startPt;
            coordsMT = coordCurr; 
            orientBank = [orientationOld];
%             totalImprovedCones = 2;
            currImproveCones = 0;
            % iteratively propagate along microtubule
            while success && iter< max_iter
                
                if iter <= 2, fov = 1.5*fieldOfVision; else, fov = fieldOfVision; end

                [coordNext, success, orientationNext] = Methods.estimateNextCurveCoord( coordCurr, orientationOld, stepSize, visibility, fov, imageMT, bkg_thresholding, plotflags);
                if ~success % try again for a limited number of times with smaller steps
                    [coordNext, success, orientationNext] = Methods.estimateNextCurveCoord( coordCurr, orientationOld, stepSize, visibility/2, fov, imageMT, bkg_thresholding, plotflags);
%                     currImproveCones = currImproveCones+1;
                end

                if success
                    orientBank = [ orientBank, orientationNext];
                    try, orientationOld = mean( orientBank(end-2:end) ); end, try, orientationOld = mean( orientBank(end-1:end) ); end
                    coordCurr = coordNext;
                    coordsMT = [ coordsMT , coordNext];

                end

                iter = iter+1;
                
            end
            
        end
        % }}}
        
        % estimateNextCurveCoord{{{
        function [coordNew,success,phiFinal] = estimateNextCurveCoord( coordOld, orientation, stepSize, visibility, fieldOfVision, helperImage, threshold_signal_based_on_environment, plotflags)

            % Given a starting point and some other parameters, this
            % function will find the next point that is connected by high
            % intensity.

            if nargin < 8
                plotSuccess = 0;
                plotFail = 0;
            else
                plotSuccess = plotflags.success;
                plotFail = plotflags.fail;
            end

            if nargin < 7, threshold_signal_based_on_environment = 1; end

            % This is the start point for propagation
            x0 = coordOld(1);
            y0 = coordOld(2);
            
            % We will create a finer image to allow us to use smaller step
            % sizes.
            [x, y] = meshgrid( 1 : size( helperImage, 1) );
            xVecNew = 1 : 0.25: size( helperImage, 1);
            [xi, yi] = meshgrid( xVecNew );
            try
                imSub = interp2(x,y,helperImage,xi,yi,'linear');
            catch
                stoph = 1;
            end
            clear x y xi yi
           
            % radial integration {{{
            angRange = deg2rad( fieldOfVision/2 );

            % Define the step sizes to use for the angular sweep in phi. We
            % start at 90 degrees to the proposed orientation of the tube
            PhiStep = deg2rad ( 5 ); % 2.5 degrees
            PhiVec = orientation - angRange : PhiStep : orientation + angRange;

            % Pre-allocate array for storing intensity data with varying phi.
            IntPhi = zeros( 1, length( PhiVec) );

            rp = visibility; % used for radial integration (allows for better angular determination)
            rmin = 1;
%             figure;
            for jPhi = 1 : length( PhiVec)
                % Extract the relevant angles for drawing a non-zero angle
                phi1 = PhiVec( jPhi) - PhiStep;
                phi2 = PhiVec( jPhi) + PhiStep;

                % Find the coordinates of the extreme points on the arc length
                % of this pizza slice.
                X1 = x0 + ( visibility * cos( [ phi1, phi2]) );
                Y1 = y0 + ( visibility * sin( [ phi1, phi2]) ); % -ve because matlab orientation is reversed( origin is at top left)

                % find minimum points at distance rmin away
                X2 = x0 + ( rmin * cos( [ phi1, phi2]) );
                Y2 = y0 + ( rmin * sin( [ phi1, phi2]) ); 

                xx = [ X1, X2];
                yy = [ Y1, Y2];

                xRd = []; yRd = [];
                for jPt = 1 : length(xx)
                    [~, xidx] = min( abs( xVecNew - round( xx(jPt), 1) ) );
                    [~, yidx] = min( abs( xVecNew - round( yy(jPt), 1) ) );
                    xRd = [ xRd, xidx ];
                    yRd = [ yRd, yidx ];
                end

                % Initialize the mask. We'll turn on the (X,Y) pixels and draw a convex
                % hull around it.
                imMask = zeros( size( imSub) );
                imMask( sub2ind( size( imSub), yRd, xRd ) ) = 1;

                % Draw a convex hull and uncover the pixel values of interest
                imMask = bwconvhull( imMask);

                imMasked = imMask .* imSub;
%                 imshowpair(imSub, imMask); pause(0.25)
                
                % Sum up values and store them
                IntPhi( jPhi) = sum( imMasked(:) ) / sum(imMask(:));
                
            end
            % We smooth the angular intensity to enable easy peak finding
            imSmooth = imgaussfilt( IntPhi, 3);
            % }}}
           
            if min(imSmooth) == max(imSmooth)
                success = 0; phiFinal = NaN; coordNew = [NaN; NaN];
%                 disp('Min Value equal to Max value in radial integration')
                plotInfo( imSmooth, [x0, y0], PhiVec, helperImage, visibility, stepSize, phiFinal, plotFail)
                return
            end
           
            % Signal Modification based on environment{{{
            imSmooth = ( max(IntPhi)-min(IntPhi) ) * mat2gray( imSmooth ) + min( IntPhi);
            
            % Now we need to come up with a threshold for defining what is
            % and what isn't a peak. We'll do local intensity search for
            % this
            if threshold_signal_based_on_environment
            
                imMask = logical(0*helperImage);
                for jY = 1 : size(helperImage, 1)
                    for jX = 1 : size(helperImage, 2)
                        if ( norm( [x0, y0] - [jX, jY] ) < visibility ) && ( norm( [x0, y0] - [jX, jY] ) > rmin )
                            imMask( jY, jX) = 1;
                        end
                    end
                end
                imBkg = helperImage(imMask);
                T = multithresh( imBkg, 1 );
                imSmooth( imSmooth < mean(imBkg) ) = mean(imBkg);
                minHeight = mean(imBkg); 
            
            else
                minHeight = mean( imSmooth);
            end
            % }}} 
            
            if min(imSmooth) == max(imSmooth)
                success = 0; phiFinal = NaN; coordNew = [NaN, NaN];
%                 disp('Min Value equal to Max value in radial integration')
                plotInfo( imSmooth, [x0, y0], PhiVec, helperImage, visibility, stepSize, phiFinal, plotFail)
                return
            end
            
            % Find the brightest peak {{{ 
            % Find Peaks Properties
            halfwidth = asin( 1/ stepSize);
            if ~isreal(halfwidth); halfwidth=asin(1.5/2); end;
            minProm = 0.0001;
            props = {'SortStr', 'descend', 'MinPeakProminence', minProm, 'MinPeakWidth', halfwidth, 'MinPeakHeight', minHeight};

            % Find the angle corresponding to the peaks in this region
            try
                [pkInt, phi0Loc, ~, ~] = findpeaks( imSmooth, props{:} );
            catch
                success = 0; phiFinal=NaN; coordNew = [NaN; NaN];
                plotInfo( imSmooth, [x0, y0], PhiVec, helperImage, visibility, stepSize, phiFinal, plotFail)
                return
            end
            phi0 = PhiVec( phi0Loc);
            
            % pick the brightest peak.
            [~, maxIdx] = max(pkInt);
            phiFinal = phi0( maxIdx);
            if length(phiFinal) > 1
                stoph = 1;
            end
            
            if ~isempty(phiFinal)
                phiFinal = phiFinal';
                x1 = x0 + ( stepSize * cos( phiFinal ) );
                y1 = y0 + ( stepSize * sin( phiFinal ) );
            else
                x1 = NaN; y1 = NaN;
            end
            
            % ensure this point is connected via intensity (disallow jumps)
            if ~isnan(x1) && ~isnan(y1)
                imJumps = 0*helperImage;
                imJumps( round(y1), round(x1) ) = 1; imJumps( round(y0), round(x0) ) = 1; imJumps = bwconvhull( imJumps);
                
                % if (x1, y1) has intensity, then there cannot be more than 1 zero pixels in the conv hull
                if helperImage( round(y1), round(x1) ) ~= 0 && sum(imJumps(:))-sum( imJumps(:).*helperImage(:) ~= 0 ) > ceil(stepSize/3)
                    x1 = NaN; y1 = NaN;
                end
                % if (x1, y1) has no intensity, then atleast 50% of the hull better have intensity.
                if sum( imJumps(:).*helperImage(:) ~= 0) < 0.5*sum( imJumps(:) )
                    x1 = NaN; y1 = NaN;
                end
            end 

            if ~isnan( x1) && ~isnan(y1)
                coordNew = [x1; y1]; 
                success = 1;
            else
                coordNew = [NaN, NaN];
                success = 0;
            end
            % }}} 
            
            plotInfo( imSmooth, [x0, y0], PhiVec, helperImage, visibility, stepSize, phiFinal, plotSuccess)

            % plotInfo {{{
            function plotInfo( angInt, coord, phiInfo, sourceImage, vis, step, phiFinal, plotflag)
                if plotflag==0
                    return
                end

                edgeX = coord(1) + (vis * cos( [phiInfo(1), phiInfo(end)] ) );
                edgeY = coord(2) + (vis * sin( [phiInfo(1), phiInfo(end)] ) );
                newX = coord(1) + (step * cos( phiFinal) );
                newY = coord(2) + (step * sin( phiFinal) );

                x0 = coord(1); y0 = coord(2); 

                figure;
                subplot(121)
                imshow( sourceImage, []); hold on;
                plot([edgeX(1), x0, edgeX(2)], [edgeY(1), y0, edgeY(2)], 'r-', 'LineWidth', 2); 
                if ~isnan( newX), plot( newX, newY, 'r*', 'MarkerSize', 10); end, hold off;

                subplot(122)
                plot(phiInfo, angInt, 'r', 'LineWidth', 3); hold on
                if ~isnan( phiFinal), plot( [phiFinal, phiFinal], [min(angInt), max(angInt)], 'r-.', 'LineWidth', 2); end, hold off;

            end
            % }}}


        end
        % }}}

        % RemoveCurveFromImage {{{
        function imPost = RemoveCurveFromImage( imPrior, coord, imMask)
            
            imPost = imPrior;
            imPost( imPost < 0) = 0;

            % Find pixel coordinates of curve
            imMT = logical( 0*imPrior);
            [cdX, cdY] = Methods.InterpolateCoords( coord(1,:), coord(2,:), 10);
            idxMT = sub2ind( size(imPrior), round(cdY), round( cdX) );
            imMT( idxMT) = 1;
            imMT = imdilate( imMT, strel('disk', 2) );
            idxMT = find( imMT);

            % Set microtubule pixels to zero
            imPost( idxMT) = 0;
            minV = min( imPost(:) ); maxV = max( imPost(:) ); 
            imPost = mat2gray( imPost) * ( maxV-minV) + minV;
            imPost = imPost .* imMask;
            imPostLog = logical(imPost); 
            imPost = bwareafilt( imPostLog, [15, Inf]) .* imPrior;
            imPost = imPost.*imMask;

        end
        %  }}}

        % RemoveSeedFromImage {{{
        function imPost = RemoveSeedFromImage( imPrior, coord, width, imMask)
            
            imPost = imPrior;
            imPost( imPost < 0) = 0;

            % de-intensity the mtoc point
            for jX = 1 : size(imPost, 2); for jY = 1 : size(imPost, 1)
                if norm( [jX, jY] - coord ) <= width 
                    imPost(jY, jX) = 0;
                end
            end; end
            imPost = imPost.*imMask;

        end
        %  }}}
        
        % InterpolateCoords {{{
        function [coordsNewX, coordsNewY] = InterpolateCoords( coordX, coordY, numIntervals)

            if any( size(coordX) ~= size(coordY) ), error( 'input size of coords is different'), end

            nC = length( coordX);

            coordsNewX = [];
            coordsNewY = [];

            for jC = 1 : nC-1
                cX = linspace( coordX(jC), coordX(jC+1), numIntervals);
                cY = linspace( coordY(jC), coordY(jC+1), numIntervals);
                coordsNewX = [ coordsNewX, cX(1:end-1)];
                coordsNewY = [ coordsNewY, cY(1:end-1)];
                if jC == nC-1
                    coordsNewX = [coordsNewX, cX(end)]; coordsNewY = [coordsNewY, cY(end)];
                end

            end

            if isempty(coordsNewX)
                coordsNewX = coordX;
                coordsNewY = coordY;
            end

        end
        % }}}
        
        % FilterImageForCurveDetection {{{
        function [imFilt, imSteer] = FilterImageForCurveDetection( imageIn, varargin)

            % Parse Arguments
            opts = parseArgs( imageIn, varargin{:});

            % Dimension of image
            dim = length(size(imageIn));

            % Apply mask to image
            imMasked = imageIn .* opts.Mask;

            % Set masked region values to median of unmasked image. 
            % This helps avoid finding lines on the edges of the mask.
            imFull = imMasked;
            imFull( opts.Mask == 0) = median( imFull( imFull ~= 0) );

            % Filtered image
            % Bandpass filtering to find features of specific width
            imFilt = 0*imageIn;
            for jZ = 1 : size( imageIn, 3)
                imFilt( :,:,jZ) = ImageData.FilterGaussBandpass( imFull(:,:,jZ), opts.WidthMax, opts.WidthMin).* opts.Mask(:,:,jZ);
            end

            % Steerable image with detected curves
            imSteer = 0*imageIn;
            for jZ = 1 : size(imageIn, 3)
                [~, ~, imSteer(:,:,jZ)] = steerableDetector( imFull(:,:,jZ), opts.WidthEstimate, sqrt(2)*1.3 );
            end

            imSteer = imSteer .* opts.Mask;

            % parseArgs {{{
            function opts = parseArgs( image, varargin)
                % Possible Arguments
                % Mask
                % WidthMin
                % WidthMax
                % WidthEstimate
                % Plot

                % Dimension of image
                dim = length( size( image) );

                % Default Values
                defaultMask = image ~= 0;
                defaultWidthMin = 1;
                defaultWidthMax = 4;
                defaultWidthEstimate = 4;
                defaultPlot = 0; % 0: no plots, 1: yes plots 

                % Input Parser
                p = inputParser;

                % Mask 
                validMask = @(x) length(size(x)) == dim && all( size(x) == size(image));
                addParameter( p, 'Mask', defaultMask, validMask);

                % Function for positive integer
                validPosInt = @(x) length(x)==1 && isnumeric(x) && x>0; 

                % WidthMin 
                addParameter( p, 'WidthMin', defaultWidthMin, validPosInt);
                % WidthMax 
                addParameter( p, 'WidthMax', defaultWidthMax, validPosInt);
                % WidthEstimate
                addParameter( p, 'WidthEstimate', defaultWidthEstimate, validPosInt);
                % Plot 
                validPlot = @(x) length(x)==1 && isnumeric(x) && (x==0 || x==1); 
                addParameter( p, 'Plot', defaultPlot, validPlot);

                parse( p, varargin{:} );
                opts = p.Results;

            end
            % }}}

        end
        % }}}

        % GetImageStats {{{
        function stats = GetImageStats( image, verbose)

            if nargin == 1
                verbose = 1;
            end

            % Dimension
            stats.Dim = length( size( image));

            % Mean, Median and Sigma
            imVals = image( image ~= 0); imVals = imVals(:);
            stats.Median = median( imVals );
            stats.Mean = mean( imVals);
            stats.Sigma = std( imVals);
            
            % Otsu threshold
            otsu = multithresh( imVals, 2); 
            stats.Otsu = otsu(2);

            % Common thresholds for features
            stats.ThreshHigh = stats.Mean + 3*stats.Sigma;
            stats.ThreshLow = stats.Mean + stats.Sigma;
            
            % Max
            stats.Max = max(image(:));
            stats.ThreshHigh = (stats.Max - stats.Median)/2 + stats.Median;
            stats.ThreshLow = (stats.Max - stats.Median)/4 + stats.Median;

            % Display stats
            if verbose 
                disp('Image Stats:')
                fprintf( 'Dimension = %.1f\n', stats.Dim);
                fprintf( 'Mean = %.2f\n', stats.Mean);
                fprintf( 'Median = %.2f\n', stats.Median);
                fprintf( 'Sigma = %.2f\n', stats.Sigma);
                fprintf( 'Otsu = %.2f\n', stats.Otsu);
                fprintf( 'Int Thresh High = %.2f\n', stats.ThreshHigh);
                fprintf( 'Int Thresh Low = %.2f\n', stats.ThreshLow);
            end

        end
        % }}}
        
        % Opts2VarsExcept {{{
        function vars = Opts2VarsExcept( opts, varargin)
            
            for jv = 1 : length( varargin)
                try
                    opts = rmfield( opts, varargin{jv});
                end
            end
            vars = Methods.struct2cellvars( opts);
            
        end
        % }}}
        
        % struct2cellvars {{{
        function vars = struct2cellvars( st)
            
            fields = fieldnames(st);
            vals = struct2cell(st);
            n = length( fields);
            
            vars = cell(1, 2*n);
            
            vars(1:2:end-1) = fields;
            vars(2:2:end) = vals;
            
        end
        % }}}
        
        % CheckEscapeMask {{{
        function [status,amt] = CheckEscapeMask( imageIn, mask, sens)
            % Check if features outside 2D mask
            % imageIn: simulated image of any number of features

            if nargin < 3
                sens = 0.2;
            end
            status = 0; amt=0;
            % Maximum Intensity of simulated image
            maxSim = max( imageIn(:) );
            
            % Intensity image outside mask
            intOutside = imcomplement( logical( mask) ) .* imageIn;

            % Intensity threshold for penalizing
            intThresh = sens * maxSim;

            if max( intOutside(:) ) > intThresh 
                status = 1;
                amt = sum( intOutside(:) );
            end
        end
        % }}}

    end
end


