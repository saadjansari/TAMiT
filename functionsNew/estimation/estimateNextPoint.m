function [coordNew,success,phiFinal] = EstimateNextPoint( coordOld, orientation, stepSize, visibility, fieldOfVision, helperImage, threshold_signal_based_on_environment, plotflags)

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
            imSub = interp2(x,y,helperImage,xi,yi,'linear');

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

                X = [ X1, X2];
                Y = [ Y1, Y2];

                xRd = []; yRd = [];
                for jPt = 1 : length(X)
                    [~, xidx] = min( abs( xVecNew - round( X(jPt), 1) ) );
                    [~, yidx] = min( abs( xVecNew - round( Y(jPt), 1) ) );
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
                if helperImage( round(y1), round(x1) ) ~= 0 && sum(imJumps(:))-sum( imJumps(:).*helperImage(:) ~= 0 ) > 2
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
