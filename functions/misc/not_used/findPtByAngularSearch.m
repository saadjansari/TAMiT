function [phiFinal, x1, y1] = findPtByAngularSearch( imPlane,startPt, r0, visibilityRatio, plotFlag, angleStart, angleRange)
% Given an 2D-image and a starting point, this function scans around in 360
% degrees and looks for peaks. It radially integrates up to some radius.
% Usage: given a point spoecified on a line/curve, this funciton can be
% used to jump to another point on the connected intensity. It can be used
% to create a bead-like model describing the line/curve

fineSearch = 0;
if nargin == 7
    fineSearch = 1;
elseif nargin == 6
    error('7th argument required')
elseif nargin == 4
    plotFlag = 0;
end
    
    

keepBest2 = 1;

% Extract initial points
x0 = startPt(1);
y0 = startPt(2);

% apply conv filter to image to re4duce noise
imC = imgaussfilt( imPlane, 1);


xVec = 1 : size( imPlane, 1);
xVecNew = 1 : 0.25: size( imPlane, 1);
[x, y] = meshgrid( xVec );
[xi, yi] = meshgrid( xVecNew );
imPlaneNew = interp2(x,y, imPlane,xi,yi,'cubic');
% imCNew = imgaussfilt( imPlaneNew, 1);
imCNew = imPlaneNew;



% Define the step sizes to use for the angular sweep in phi
PhiStep = deg2rad ( 2.5 ); % 1 degrees
PhiVec = - pi : PhiStep : pi-PhiStep;

% Pre-allocate array for storing intensity data with varying phi.
IntPhi = zeros( 1, length( PhiVec) );


% figure; set(gcf, 'pos', get(0, 'ScreenSize') );
rp = visibilityRatio *r0; % used for radial integration (allows for better angular determination
rmin = 2; 


% figure; hold on;
for jPhi = 1 : length( PhiVec)
    % Extract the relevant angles for drawing a non-zero angle
    phi1 = PhiVec( jPhi) - PhiStep/2;
    phi2 = PhiVec( jPhi) + PhiStep/2;

    % Find the coordinates of the extreme points on the arc length
    % of this pizza slice.
    X1 = x0 + ( rp * cos( [ phi1, phi2]) );
    Y1 = y0 + ( rp * sin( [ phi1, phi2]) );

    % find intermediate points
    X2 = x0 + ( rmin * cos( [ phi1, phi2]) );
    Y2 = y0 + ( rmin * sin( [ phi1, phi2]) );
    
    % Append the SPB to get the third point for the convex hull
%     X = [ x0, X1];
%     Y = [ y0, Y1];
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
    imMask = zeros( size( imCNew) );
    imMask( sub2ind( size( imCNew), yRd, xRd ) ) = 1;

    % Draw a convex hull and uncover the pixel values of interest
    imMask = bwconvhull( imMask);

    imMasked = imMask .* imCNew;
%     imagesc(imMask); axis equal; pause(0.025)    

    % Sum up values and store them
    IntPhi( jPhi) = sum( imMasked(:) ) / sum(imMask(:));
    
end
% hold off


% Smooth this rotational intensity and look for a global minmum, we
% will use this minimum to define a new vector of phi values (for
% continuity in peaks)
imSignal = mat2gray( imgaussfilt( IntPhi, 3) );
[~, PhiMinIdx] = findpeaks( -imSignal, 'NPeaks', 1, 'SortStr', 'descend');
% offset the intensity pattern and phi vector so the minimum is on the
% endpoints
PhiVecNew = [ PhiVec( PhiMinIdx : end), PhiVec(1: PhiMinIdx-1)+2*pi ];
IntPhiNew = [ IntPhi( PhiMinIdx : end), IntPhi(1: PhiMinIdx-1) ];
imSignalNew = [ imSignal( PhiMinIdx : end), imSignal(1: PhiMinIdx-1) ];

imSignalC = ( max(IntPhiNew)-min(IntPhiNew) ) * mat2gray( imSignalNew ) + min( IntPhiNew);

% We need a threshold intensity value, and instead of finding a threshold
% for the whole image, we'll just use the threshold from a local region
imMask = 0*imPlane;
for jY = 1 : size(imPlane, 1)
    for jX = 1 : size(imPlane, 2)
        if ( norm( startPt - [jX, jY] ) < rp ) && ( norm( startPt - [jX, jY] ) > rmin )
            imMask( jY, jX) = 1;
        end
    end
end
% figure; imagesc(imMask); axis equal
imMasked = imMask .* imPlane;
imVals = imMasked( imMasked > 0);
T = multithresh( imVals, 1 );
imSignalC( imSignalC < median(imVals) ) = median(imVals);
if T < median(imVals); T = median(imVals); end

% Find Peaks Properties
halfwidth = asin( 1.5/ r0);
minProm = 0.01;
minHeight = mean( [T, median(imVals)]);
props = {'SortStr', 'descend', 'MinPeakProminence', minProm, 'MinPeakWidth', halfwidth, 'MinPeakHeight', median(imVals)};

% Find the angle corresponding to the peaks in this region
% [~, phi0, w, ~] = findpeaks( mat2gray(imSignal), PhiVec, 'SortStr', 'descend', 'MinPeakWidth', halfwidth, 'MinProminence', minProm);
[pkInt, phi0Loc, w, p] = findpeaks( imSignalC, props{:} );
phi0 = PhiVecNew( phi0Loc);

if fineSearch
    angLim = angleStart + [1,-1]*angleRange;
else
    angLim = [pi -pi];
end

% re-map phi to lie between - pi to +pi
phi0 = phi0 - (phi0 > pi) *2*pi;

% Now we'll begin the peak selection procedure. Is there a expected range
% for a peaK?
% if any(angLim < -pi)
%     % map angles to -2pi to 0
%     phi0 = phi0 - pi;
%     phiFinal = phi0( (phi0 < angLim(1) ) & ( phi0 > angLim(2) ) ) + pi;
% elseif any(angLim > pi)
%     % map angles to 0 to +2pi
%     phi0 = phi0 + pi;
%     phiFinal = phi0( (phi0 < angLim(1) ) & ( phi0 > angLim(2) ) ) - pi;
% else
%     phiFinal = phi0( (phi0 < angLim(1) ) & ( phi0 > angLim(2) ) );
% end
% 
% if keepBest2
%     phiFinal(3:end) = [];
% end


% Check if the angles in phi0 lie between the angle Limits (angles are
% periodic so something bigger is needed here
% Angle limits can exceed -pi and +pi

if angLim(2) < -pi
    
    % check if phi0 lies between -pi and upper angle limit
    log1 = phi0 > -pi & phi0 < angLim(1);
    
    % check if phi0 lies between angle lower limit and +pi (since angles
    % are modded)
    log2 = phi0 > angLim(2)+2*pi & phi0 < pi;
    
    phiFinal = phi0( log1 | log2);
    pkIntFinal = pkInt( log1 | log2);
    
elseif angLim(1) > pi
    
    % check if phi0 lies between -pi and upper angle limit
    log1 = phi0 > -pi & phi0 < angLim(1)-2*pi;
    
    % check if phi0 lies between angle lower limit and +pi (since angles
    % are modded)
    log2 = phi0 > angLim(2) & phi0 < pi;
    
    phiFinal = phi0( log1 | log2);
    pkIntFinal = pkInt( log1 | log2);
    
else
    phiFinal = phi0( (phi0 < angLim(1) ) & ( phi0 > angLim(2) ) );
    pkIntFinal = pkInt(  (phi0 < angLim(1) ) & ( phi0 > angLim(2) ) );
end

% if more than 1 peak while using fineSearch, only pick the brightest one.
if fineSearch
    [~, maxIdx] = max(pkIntFinal);
    phiFinal = phiFinal( maxIdx);
end

% Find the corresponding point in cartesian space using r0 and the peak
% angle
if ~ isempty(phiFinal)
    phiFinal = phiFinal';
    x1 = x0 + ( r0 * cos( phiFinal ) );
    y1 = y0 + ( r0 * sin( phiFinal ) );

    ptNext = [x1 , y1];

    % plot image showing starting point, over circle, and best link points

    % create boundary for circle
    th = 0:pi/50:2*pi;
    xCir = r0 * cos(th) + x0;
    yCir = r0 * sin(th) + y0;

    if plotFlag
        
        h = figure; 
        subplot(131);
        findpeaks( imSignalC, PhiVecNew, props{:}, 'Annotate', 'extents' );
        if fineSearch
            hold on; plot( [angLim(2), angLim(2)], [min(imSignalC), max(imSignalC)], 'r:')
            plot( [angLim(1), angLim(1)], [min(imSignalC), max(imSignalC)], 'c:')
        end
        hold off
        
        
        % create figure 1 of local angular search
        % display everything
        msz = 10; lwsz = 2;
        subplot(132)
        imagesc( imPlane); colormap gray; axis equal; set(gca, 'xlim', [1 size(imPlane,1)], 'ylim', [1 size(imPlane,1)]); hold on; %display image
        if fineSearch
            th2 = angLim(2) : pi/100: angLim(1);
            xCir2 = r0 * cos(th2) + x0;
            yCir2 = r0 * sin(th2) + y0;
            plot(xCir, yCir, 'r','LineWidth', lwsz); % plot circle boundary
            plot(xCir2, yCir2, 'g','LineWidth', lwsz); % plot circle boundary
        else
            plot(xCir, yCir, 'g','LineWidth', lwsz); % plot circle boundary
        end
        plot( x0, y0, 'bx', 'MarkerSize', msz, 'LineWidth', lwsz); % plot start point
        for jPt = 1 : size( ptNext, 1)
            plot( ptNext( jPt, 1), ptNext( jPt, 2), 'bo', 'MarkerSize', msz, 'LineWidth', lwsz); % plot possible next links
        end
        % create legend entry
        LH(1) = plot(nan, nan, 'bx', 'MarkerSize', msz); L{1} = 'Start Point';
        LH(2) = plot(nan, nan, 'bo', 'MarkerSize', msz); L{2} = 'Next Point';
        LH(3) = plot(nan, nan, 'r-', 'MarkerSize', msz); L{3} = 'Forbidden Angles';
        LH(4) = plot(nan, nan, 'g-', 'MarkerSize', msz); L{4} = 'Allowed Angles';
        legend(LH, L);
        set(gca, 'FontSize', 15);
        title('Result of Angular Search')
        hold off

        
        % create Figure 2 showing connection between links
        msz = 6; lwsz = 2;
        subplot(133);
        imagesc( imPlane); colormap gray; axis equal; hold on; set(gca, 'xlim', [1 size(imPlane,1)], 'ylim', [1 size(imPlane,1)]); %display image
        plot( [ x0, ptNext(:,1)'], [ y0, ptNext(:,2)'], 'bx', 'MarkerSize', msz, 'LineWidth', lwsz); % plot points
        for jPt = 1 : size( ptNext, 1)
            plot( [ x0, ptNext(jPt,1) ], [ y0, ptNext(jPt,2)], 'b-', 'LineWidth', lwsz); % plot line connecting points
        end
        hold off
        set(gcf, 'WindowState', 'maximized');
    end
    
else
    x1 = NaN; y1 = NaN; phiFinal = 0;
end

end

