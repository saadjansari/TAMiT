    function boneProps = analyzeLineImage( bone)
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
