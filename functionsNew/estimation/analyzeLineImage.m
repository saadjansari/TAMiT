    function lineProps = analyzeLineImage( lineImage, lineCoords)
        % findCurveVariance: find the variance of the local curvatures of
        % an object as well as some other properties
        %  if providing a lineImage, use syntax analyzeLineImage( lineImage, [])
        %  if providing coordinates, use syntax analyzeLineImage( lineImage, lineCoords) note in this case lineImage is only usedwith size() function

        if nargin ==1 || isempty(lineCoords)

            % find endpoint of line
            Line = logical( lineImage);
            bw = bwmorph( Line, 'endpoints');
            endPtIdx = find(bw == 1);
            [ y1, x1] = ind2sub( size(Line), endPtIdx);
            
            % pick an endpoint and measure geodesic distance to each other
            % pixel to get order list of connected pixels
            mat_dist = bwdistgeodesic( Line, x1(1), y1(1), 'quasi-euclidean'); %'quasi-euclidean' for 8-connectivity
            comp = find( Line);
            comp(:,2) = mat_dist(comp(:,1));
            ordered_list_ind = sortrows(comp,2);
        
        else
try
            ordered_list_ind = sub2ind( size(lineImage), lineCoords(2,:), lineCoords(1,:) ); 
catch
size(lineImage)
size(lineCoords)
end
        end

        [xyList(:,1), xyList(:,2)] = ind2sub( size(lineImage), ordered_list_ind );     
        numPix = size( xyList,1);

        % We will measure the following properties:
        %  1) global orientation
        %  2) local orientations
        %  3) length
        %  4)
        %  5)
        %  6)

        % Find Local Orientations: use a moving filter of 5 pixels
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
            
            thetaNextPossible = mean( atan2(yD, xD) ) + [-2*pi, 0, +2*pi];
            if ~isempty( thetaList)
                [~, thetaNextIdx] = min( abs( thetaNextPossible - thetaList(end) ) );
            else
                thetaNextIdx = 2;
            end
            thetaList = [ thetaList, thetaNextPossible( thetaNextIdx) ];
            
        end
        
        % Find Global Orientation
        thetaGlobal = mean( thetaList);

        % Find Path Length
        len = numPix;

        % store information about line/curve
        lineProps.length = len;
        lineProps.theta = thetaGlobal;
        lineProps.thetaVar = std( thetaList);
        lineProps.thetaList = thetaList;
        lineProps.orderedPixelList = xyList;
        lineProps.orderedIdxList = ordered_list_ind;
end
