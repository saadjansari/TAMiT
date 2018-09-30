function finalNetwork = reconstructNetwork( brokenNetwork, intensityImage, params)
% reconstructNetwork: Take an image of a network of lines/curves that has
% some broken links. This function looks at each endpoint and finds the
% cost of connecting it to a nearby object based on local curvature
% matching, connection distance and line-axis shift matching.
% Endpoints-Endpoint connections are encouraged where possible.

plotFlag = 0;
% We start with an image of a broken network.
% We have the barebones of the network, but it needs to be reinstated.
% For each object we will determine if it should be linked to any of its
% neighbors.
disp('Reinstating the network...')

% count the number of objects;
cc = bwconncomp( brokenNetwork); nObj = cc.NumObjects;
fprintf('Total Number of objects = %d\n', nObj) 

% The first step is to determine which pixels are end-pixels. This is based
% on adding up values in a 3x3 neighborhood and looking at distances of
% positive values to conclude if a pixel is a line-center pixel or an
% end-pixel
endPoints = findNetworkEndpoints( brokenNetwork);

% display image of endpoints
if plotFlag
    imEnd = 0 * brokenNetwork;
    imEnd( [endPoints.pixId] ) = 1;
    imEnd2 = zeros( [ size(imEnd), 3] );
    imEnd2(:,:,1) = brokenNetwork-imEnd;
    imEnd2(:,:,2) = imEnd;
    dispImg( imEnd2); title('Endpoints of Network in Green'); set(gca, 'FontSize', 15);
end

% Now we have the endpoints of all the objects. We'll find the local
% orientation of all these endpoints (measured from the positive
% horizontal axis to the positive axis anticlockwise).
endPoints = findEndpointOrientation( endPoints);

% Now, we have the local orientation of all the endpoints.
% Next, we will draw a cone from each endpoint in its local orientation and
% find objects it can be linked with. We will also record information about
% the objects it can be linked with
linkPoints = findLinkPoints( endPoints, params.costs);

% Now, we will find the cost of every possible link
linkCosts = findLinkCosts( linkPoints, intensityImage, params.costs);

linkCosts = removeDuplicateLinks( linkCosts);

% establish links whose cost is below maximum cost
finalNetwork = completeLinks( brokenNetwork, linkCosts);

% display image of broken and final network
if plotFlag
    imDisp = zeros( [size(brokenNetwork), 3] );
    imDisp( :,:,1) = brokenNetwork; imDisp( :, :, 2) = finalNetwork-brokenNetwork;
    dispImg( imDisp); title('Reinstated Network: Greens are new connections'); set(gca, 'FontSize', 15);
end

disp('Network successfully reinstated')
%% endPoints = findNetworkEndpoints( networkIn) {{{
    function endPoints = findNetworkEndpoints( networkIn)
        % findNetworkEndpoints: finds the endpoints of a line/curve or any
        % shape composed of lines/curves. endPoints output is an mx4 array
        % where m is the number or endpoints. The first column gives the
        % endpoint pixel index, 2nd and 3rd column give the x(horz) and y(vert)
        % pixel coordinate. The 4th column gives the object number that
        % this endpoint belongs to. endPoints is now a strucutre with
        % similar layout but fields 'pixId, x, y, objId'.
        
        % for each object, lets find the end pixel. The end pixel is such that if
        % we add up a 3x3 neighborhood, we get a value of 2. We'll first try a
        % connectivity of 4 and see if a pixel has a neighborhood sum of 2. This
        % will indicate a vertical/horizontal/L connection. If there isnt such a
        % connection, we'll look at the diagonal elements and see if the sum of
        % pixels is 1. ! indicates that there is just a single diagonal link to one
        % side of the object.
        endPoints = [];
        
%         % for each object in the broken network
%         for jObj = 1 : nObj
%             
%             % get all its pixels, and turn on just the object
%             pixList = cc.PixelIdxList{ jObj};
%             imObj = 0 * brokenNetwork;
%             imObj( pixList) = 1;
%             % get cartesian indices of pixel
%             [yList, xList] = ind2sub( size(imObj), pixList );
%             
%             % check if any pixel satisfies our criterion for being an
%             % endpixel.
%             for jPix = 1 : length(pixList)
%                 
%                 y = yList( jPix); x = xList( jPix);
%                 
%                 % get the local 3x3 neighborhood
%                 nhood = imObj( y-1:y+1, x-1:x+1);
%                 
%                 % apply criterion for endpoint
%                 if sum( nhood( 2:2:end) ) < 2
%                     % if there aren't any horizontal/vertical or L-shaped
%                     % connections, it is an endpoint if there is only one
%                     % turned on pixel in the neighborhood. It is also an
%                     % endpoint if there are 2 turned on pixels but they are
%                     % very close together, indicating that something ends
%                     % at the current pixel.
%                     
%                     if sum( nhood(:) ) - 1 == 1
%                         
%                         endPoints = [ endPoints; pixList( jPix), x, y, jObj];
%                         
%                     elseif sum( nhood(:) ) - 1 == 2
% 
%                         % find the pixels that are turned on
%                         idx = find( nhood);
%                         idx = idx( idx ~= 5); % not the current pixel
%                         [yi, xi] = ind2sub( size(nhood) , idx );
% 
%                         % find dist between pixels
%                         if norm([diff(xi) diff(yi)]) < 2.1
%                             
%                             endPoints = [ endPoints; pixList( jPix), x, y, jObj];
%                             
%                         end
%                         
%                     end
%                     
%                 end
%                 
%             end
% 
%         end
%         
%         for jEnd = 1 : size( endPoints, 1)
%             
%             endPointsStruct( jEnd).pixId = endPoints( jEnd, 1);
%             endPointsStruct( jEnd).x = endPoints( jEnd, 2);
%             endPointsStruct( jEnd).y = endPoints( jEnd, 3);
%             endPointsStruct( jEnd).objId = endPoints( jEnd, 4);
%             
%         end
%         
%         endPoints = endPointsStruct;

        imEndPoints = bwmorph( networkIn, 'endpoints');
        
        % for every connected pixels, find its centroid and that will be
        % the endpoint
        ccc = regionprops( imEndPoints, 'Centroid');
        for jEnd = 1 : length(ccc)
            
            endPointsStruct( jEnd).x = round( ccc(jEnd).Centroid(1) );
            endPointsStruct( jEnd).y = round( ccc(jEnd).Centroid(2) );
            endPointsStruct( jEnd).pixId = sub2ind(size(imEndPoints), endPointsStruct( jEnd).y, endPointsStruct( jEnd).x);
            
            % get the object number associated with this endpoint
            objs = cellfun( @(v) any( v(:) == endPointsStruct( jEnd).pixId), cc.PixelIdxList );
            endPointsStruct( jEnd).objId = find( objs);
            if isempty( find( objs) )
                error('Oh why o why wont...you fix me')
            end
        end
        
        endPoints = endPointsStruct;
        
    end
%% }}}

%% endPoints = findEndpointOrientation( endPoints) {{{
    function endPoints = findEndpointOrientation( endPoints)
        % findEndpointOrientation : finds the local orientation of
        % endPoints by looking at the closest pixels connected to the
        % endpoint.
        
        % create a solid circle of radius 5 pixels that will be centered at
        % each endpoint to find the closest pixels.
        searchRad = 7;
        imCirc = zeros( searchRad*2+1);
        centerCirc = [ searchRad+1, searchRad+1];
        for jX = 1 : size( imCirc, 2)
            for jY = 1 : size( imCirc, 1)
                if norm( centerCirc-[jY,jX] ) < searchRad
                    imCirc( jY, jX) = 1;
                end
            end
        end

        % Center the circle at each endpoint, find its closest pixels and
        % then determine the local orientation by running reigonprops

        for jEnd = 1: length( endPoints)
            
            % find the object this endpoint belongs to, and create an image
            % with just the object turned on.
            jObj = endPoints( jEnd).objId;
            imObj = 0 * brokenNetwork;
            imObj( cc.PixelIdxList{ jObj} ) = 1;
            
            % create an image with just the endpoint turned on, which
            % allows us to put on the circle here.
            imEnd = 0 * imObj;
            imEnd( endPoints(jEnd).pixId ) = 1;
            imEndPad = padarray( imEnd, [searchRad, searchRad], 0); % temporarily pad it to ensure no missing index errors
            
            % create image with circle centered at endpoint
            imEndCirc = 0 * imEndPad; 
            [y, x] = find( imEndPad );
            imEndCirc( y-searchRad : y+searchRad, x-searchRad : x+searchRad) = imCirc;

            % remove padding and generate image with endpoint's closest
            % pixels turned on.
            imEndCirc = imEndCirc( searchRad+1: end-searchRad , searchRad+1: end-searchRad );
            imEndLocal = imEndCirc .* imObj;

            % find the local orientation
            st = regionprops( imEndLocal, 'Centroid', 'Orientation');
            phi = deg2rad(st.Orientation);
            loc = st.Centroid;
            [y, x] = find( imEnd );
            
            % We need to make a coordinate shift at this point. regionprops
            % is a weird function in that it gives orientation measurements
            % from an axis that can change. Positive and negative phi value
            % are measure from positive and negative horizontal axes. We
            % will transform to a coordinate system where phi will be
            % measured from the positive x axis (starting at 0 radians) and
            % sweeping toward the positive y axis (in images, +y posints
            % down). phi runs from 0 to 2pi radians.
            if phi > 0
                phi = 2*pi - abs(phi);
            elseif phi < 0
                phi = pi + abs(phi);
            elseif phi == 0
                phi = pi;
            end
            
            % because the way phi is measured, we only ever sweep 2
            % quadrants. We will fix this because the orientation of an
            % endpoint can be either along phi or opposite phi.
            % We will draw a vector from the centroid in the direction of
            % phi, and if the vector seems to go away from the endpoint,
            % then we will reverse phi
            director = loc + 2*[ cos(phi), sin(phi)];
            if norm( director - [x,y] ) > norm( loc-[x,y] )
                phi = phi-pi;
            end
            
            % store local orientation
            endPoints( jEnd).phi = phi;
            
        end

    end
%% }}}

%% linkPoints = findLinkPoints( endPoints, cost) {{{
    function linkPoints = findLinkPoints( endPoints, cost)
        % findLinkPoints: finds the points that can be linked by drawing a
        % cone at each endpoint in the direction of its local orientation
        
        stats = regionprops( cc, 'Area');
        
        lengthCone = cost.maxDistLink;
%         phiSpreadCone = atan( cost.maxDistLinkPerp / cost.maxDistLink);
        phiSpreadCone = deg2rad(30);
        linkPoints = [];
    
        for jEnd = 1 : length( endPoints )
            
            % load values associated with current endpoint
            jObj = endPoints( jEnd).objId;
            idxEnd = endPoints( jEnd).pixId;
            xEnd = endPoints( jEnd).x;
            yEnd = endPoints( jEnd).y;
            phi = endPoints( jEnd).phi;
            
            % Create image with all objects except the one belonging to our
            % endpoint.
            imOther = brokenNetwork;
            imOther( cc.PixelIdxList{ jObj} ) = 0;

            % Create a mask with just the search cone region turned on
            phiSweep = linspace( phi-phiSpreadCone, phi+phiSpreadCone, 50);
            imCone = 0 * brokenNetwork; 
            imCone( yEnd, xEnd) = 1;
            for jPhi = phiSweep
                xV = round( round( xEnd + lengthCone * cos( jPhi) ) ); xV(xV < 1) = 1; xV( xV > size(imCone,2) ) = size(imCone,2);
                yV = round( yEnd + lengthCone * sin( jPhi)); yV(yV < 1) = 1; yV( yV > size(imCone,1) ) = size(imCone,1);
                imCone( yV, xV ) = 1;
            end
            imCone = bwconvhull( imCone);
            
            % Now, generate image containing any objects that appear in the
            % search cone.
            try
            imSearch = imCone .* imOther;
            catch
                err =1;
            end
            % We will get the pixel Indices of objects in the cone and
            % compare it with pixel indices of the brokenNetwork image to
            % determine which objects these belong to
            ccLinks = bwconncomp( imSearch);
            stLinks = regionprops( ccLinks, 'Orientation');
            numLinks = ccLinks.NumObjects;

            % for each possible link, we'll store the object it belongs to,
            % whether its an endpoint, the actual coordinate of the link
            % point, and the local orientation of the object to be linked
            % to
            
            for jLink = 1 : numLinks

                % check whether link object contains an endpoint.
                pixLink = ccLinks.PixelIdxList{ jLink};
                id = find( ismember( [endPoints.pixId], pixLink ) );

                if isempty( id)
                    % if link object does not contain an endpoint
                    
                    % specify that it does not contain an endpoint
                    endPtFlag = 0;
                    
                    % find the local orientation of the link object
                    phiLink = deg2rad( stLinks( jLink).Orientation);
                    phiDiff = min( abs( (phiLink + phi) - pi*[-4:4]) );
                    
                    % find the index of the point that should be considered
                    % for linking. This is found by extending a line based
                    % on the local orientation of the search point and
                    % finding the point closest to it.
                    imLine = 0* brokenNetwork;
                    for jR = linspace( 0, lengthCone, 2 * lengthCone)
                        imLine( round( yEnd + jR * sin( phi)), round(xEnd + jR * cos( phi) ) ) = 1;
                    end

                    % find pixel in link object with minimum perpendicular distance to imLine
                    [yLine, xLine] = find( imLine);
                    distLinkPerp = Inf;
                    linkPtX = NaN;
                    linkPtY = NaN;
                    
                    [yList, xList] = ind2sub( size( brokenNetwork), pixLink );
                    
                    % compare distance of pixels, update dist, x and y points
                    % if distance is smaller than current stored value.
                    for jLinkIdx = 1: length( pixLink)
                        for jLineIdx = 1 : length(yLine)
                            if norm( [ yList( jLinkIdx) - yLine(jLineIdx), xList( jLinkIdx) - xLine(jLineIdx) ]) < distLinkPerp
                                linkPtX = xList( jLinkIdx);
                                linkPtY = yList( jLinkIdx);
                                distLinkPerp = norm( [ yList( jLinkIdx) - yLine(jLineIdx), xList( jLinkIdx) - xLine(jLineIdx) ]);
                            end
                        end            
                    end
                    
                    % find linking distance
                    distLink = norm( [linkPtY - yEnd, linkPtX - xEnd] );
                    
                    idxLink = sub2ind( size(brokenNetwork), linkPtY, linkPtX);
                    
                    % find areas of the two objects being considered for
                    % linking
                    imLabels = bwlabel(brokenNetwork);
                    jObjLinkedTo = imLabels( idxLink);
                    areaOfLinkObject = stats( jObjLinkedTo).Area;
                    areaOfSearchObject = stats( jObj).Area;
                
                    linkPoints = [ linkPoints; idxEnd, idxLink, phiDiff, distLink, distLinkPerp, areaOfSearchObject, areaOfLinkObject, endPtFlag];

                else
                    % if link object contains a single or multiple
                    % endpoints.
                    % specify that it contains an endpoint
                    endPtFlag = 1;

                    % obtain properties of link objects since they contain
                    % endpoints
                    idxLink = [ endPoints( id).pixId];
                    objLink = [ endPoints( id).objId];
                    phiLink = [ endPoints( id).phi];
                    [ linkPtY, linkPtX ] = ind2sub( size( brokenNetwork), idxLink);
                    
                    % if more than one endpoint from the same link obejct,
                    % pick the closest one
                    
                    if length( id) > 1
                        % two endpoints of the same object
                        % pick closest one
                        d1 = norm( [yEnd - linkPtY(1), xEnd - linkPtX(1)]);
                        d2 = norm( [yEnd - linkPtY(2), xEnd - linkPtX(2)]);
                        if length(id) >= 3
                            warning('not prepared to handle 3 or more endpoints from the same object')
                        end
                        id = id(2) * (d1 > d2) + id(1) * (d1 < d2);
                    end
                    idxLink = [ endPoints( id).pixId];
                    objLink = [ endPoints( id).objId];
                    phiLink = [ endPoints( id).phi];
                    [ linkPtY, linkPtX ] = ind2sub( size( brokenNetwork), idxLink);
                    
                    % find differenc in phis to be used for cost
                    % calculation. we account for the fact that two
                    % endpoints are a perfect match if their phi values are
                    % offset by pi radians. The perfect match corresponds
                    % to a value of 0 in phiDiff below.
                    phiDiff = min( abs( (phiLink - phi) + [-pi, +pi] ) );
                    
                    % find perpendicular distance of linking the two
                    % endpoints. generate a line by extending th elocal
                    % orientation of the search point and finding the
                    % minimum distance to the linkPoint
                    imLine = 0* brokenNetwork;
                    for jR = linspace( 0, lengthCone, 2 * lengthCone)
                        imLine( round( yEnd + jR * sin( phi)), round(xEnd + jR * cos( phi) ) ) = 1;
                    end

                    % find pixel in link object with minimum perpendicular distance to imLine
                    [yLine, xLine] = find( imLine);
                    distLinkPerp = Inf;
                    % compare distance of pixels, update dist
                    % if distance is smaller than current stored value.
                    for jLineIdx = 1 : length(yLine)
                        if norm( [ linkPtY - yLine( jLineIdx), linkPtX - xLine(jLineIdx) ]) < distLinkPerp
                            distLinkPerp = norm( [ linkPtY - yLine( jLineIdx), linkPtX - xLine(jLineIdx) ]);
                        end
                    end
                    
                    % find linking distance
                    distLink = norm( [linkPtY - yEnd, linkPtX - xEnd] );
                    
                    % find areas of the two objects being considered for
                    % linking
                    imLabels = bwlabel( brokenNetwork);
                    jObjLinkedTo = imLabels( idxLink);
                    areaOfLinkObject = stats( objLink).Area;
                    areaOfSearchObject = stats( jObj).Area;
                
                    linkPoints = [ linkPoints; idxEnd, idxLink, phiDiff, distLink, distLinkPerp, areaOfSearchObject, areaOfLinkObject, endPtFlag];
                end        

            end
        
        end
           
        % create a structre to store all this data in
        for jLink = 1 : size(linkPoints, 1)
            linkPointsStruct( jLink).idx1 = linkPoints( jLink, 1);
            linkPointsStruct( jLink).idx2 = linkPoints( jLink, 2);
            linkPointsStruct( jLink).phiDiff = linkPoints( jLink, 3);
            linkPointsStruct( jLink).distLink = linkPoints( jLink, 4);
            linkPointsStruct( jLink).distLinkPerp = linkPoints( jLink, 5);
            linkPointsStruct( jLink).linkObjectSize1 = linkPoints( jLink, 6);
            linkPointsStruct( jLink).linkObjectSize2 = linkPoints( jLink, 7);
            linkPointsStruct( jLink).EndtoEndLink = linkPoints( jLink, 8);
        end
        
        if exist('linkPointsStruct')
            linkPoints = linkPointsStruct;
        end
        
    end
%% }}}

%% linkCosts = findLinkCosts( linkPoints, imPlane, costs) {{{
    function linkCosts = findLinkCosts( linkPoints, imPlane, costs)
        % findLinkCosts: finds the cost of each link
        
        linkCosts = linkPoints;
%         linkCosts = rmfield(linkCosts, 'phiDiff', 'distLink', 'diskLinkperp', 'EndtoEndLink');
        for jLink = 1 : length( linkPoints)
            
            % extract values for link to determine cost
            phi = linkPoints( jLink).phiDiff;
            dist = linkPoints( jLink).distLink;
            distPerp = linkPoints( jLink).distLinkPerp;
            ee = linkPoints( jLink).EndtoEndLink;
            area1 = linkPoints( jLink).linkObjectSize1;
            area2 = linkPoints( jLink).linkObjectSize2;
            idx1 = linkPoints( jLink).idx1;
            idx2 = linkPoints( jLink).idx2;
            
            % extract parameters that determine max values
            maxCost = costs.maxCost;
            maxDistLink = costs.maxDistLink;
            maxDistLinkPerp = costs.maxDistLinkPerp;
            maxPhiDiff_EE = costs.maxPhiDiff_EE; 
            maxPhiDiff_EL = costs.maxPhiDiff_EL;
            EEvsELScalingFactor = costs.EEvsELScalingFactor;
            linkObjectSizeScalingFactor = costs.linkObjectSizeScalingFactor;
            linkObjectSizeForScaling = costs.linkObjectSizeForScaling;
            
            
            % Phi Cost (quadratic)
            % normalization constant for phi equation will scale with dist            
            if ee == 1
                phi0 = 0.5*maxPhiDiff_EE*(1 + dist/maxDistLink);
                costPhi = ( 1/ phi0)^2 * phi^2;
            elseif ee == 0
                phi0 = 0.5*maxPhiDiff_EL*(1 + dist/maxDistLink);
                costPhi = ( 1/ phi0)^2 * phi^2;
            end
            
            
            % phi
            % Dist Cost (linear)
%             costDist = ( 1/ maxDistLink)^2 * dist^2;
            costDist = 0;
            
            % Dist perp Cost ( linear)
            costDistPerp = ( 1/ maxDistLinkPerp)^2 * distPerp^2;
            
            % add up costs and scale if necessary
            cost = costPhi + costDist + costDistPerp;
            if ee == 1; cost = cost * EEvsELScalingFactor; end
            if area1 < linkObjectSizeForScaling; cost = cost * linkObjectSizeScalingFactor; end
            if area2 < linkObjectSizeForScaling; cost = cost * linkObjectSizeScalingFactor; end
            
            % store costs
            linkCosts( jLink).cost = cost;
            linkCosts( jLink).completeLink = (cost <= maxCost);
            linkCosts( jLink).meanEndIntensity = mean( imPlane( [idx1, idx2]) );
            
            % find mean intensity of link
            imLink = logical( 0*imPlane);
            imLink( [idx1, idx2] ) = 1; imLink = bwconvhull( imLink);
            imLink = imLink .* imPlane;
            linkCosts( jLink).meanlinkIntensity = mean( imLink(imLink>0) );
            
            if linkCosts( jLink).meanlinkIntensity > linkCosts( jLink).meanEndIntensity ||...
                abs(linkCosts( jLink).meanlinkIntensity - linkCosts( jLink).meanEndIntensity)*100/linkCosts( jLink).meanEndIntensity < 5
                linkCosts( jLink).cost = linkCosts( jLink).cost *0.75;
                linkCosts( jLink).completeLink = (linkCosts( jLink).cost <= maxCost);

            end
        end
            
    end
%% }}}

%% linkCosts = removeDuplicateLinks( linkCosts) {{{
    function linkCosts = removeDuplicateLinks( linkCosts)
        
        if isempty(linkCosts)
            return
        end
        idx1 = [linkCosts.idx1];
        idx2 = [linkCosts.idx2];
        cost = [linkCosts.cost];
        
        idxRm = [];
        for jLink1 = 1 : length(idx1)
            % check if the reverse link is present in any of the other links.
            id1 = idx1( jLink1); id2 = idx2( jLink1);
            for jLink2 = 1 : length(idx1)
                if jLink2 ~= jLink1 && idx1( jLink2)==id2 && idx2( jLink2)==id1
                    c1 = cost( jLink1);
                    c2 = cost( jLink2);
                    if c1 < c2
                        idxRm = [ idxRm, jLink2];
                    elseif c1 > c2
                        idxRm = [ idxRm, jLink1];
                    elseif c1 == c2
                        idxRm = [ idxRm, jLink1*(jLink1<jLink2) + jLink2*(jLink1>jLink2)];
                    end
                end
            end
        end
        
        idxRm = unique(idxRm);
        % Now remove the duplicate links
        linkCosts( idxRm ) = [];
        
        % We also ensure that at max two links can be made to any linkpoint
        % irrespective of whether it is in index1 or index2
        idx1 = [linkCosts.idx1];
        idx2 = [linkCosts.idx2];
        cost = [linkCosts.cost];
        idxRm = [];
        numKeep = 2;
        for jLink = 1 : length(idx1)
            if isempty( idxRm) || any( jLink ~= idxRm)
                id = idx1( jLink);
                idDup = find( idx2 == id | idx1 == id);
                c = cost(idDup);
                [Csorted, sortIdx] = sort(c);
                try;idxKeep = sortIdx(1:numKeep); catch; idxKeep = sortIdx(1); end
%                 [~, idxKeep] = min( c);
                idxRm = unique([idxRm, setdiff( idDup, idDup(idxKeep) ) ]);
            end
        end
        
        idxRm = unique(idxRm);
        % Now remove the duplicate links
        linkCosts( idxRm ) = [];
        
    end
%% }}}

%% finalNetwork = completeLinks( brokenNetwork, linkCosts) {{{
    function finalNetwork = completeLinks( brokenNetwork, linkCosts)
        % completeLinks: completes the Links connecting indices in
        % linkCosts.idx1 and linkCosts.idx2 if linksCosts.completeLinks =
        % 1. the final network is a sum of the completed links and the
        % original broken network
        
        finalNetwork = brokenNetwork;
        if isfield(linkCosts, 'idx1')
            uniqueEnds = unique([linkCosts.idx1]);
        else
            uniqueEnds = [];
        end
        
        for jLink = 1 : length( uniqueEnds)
            
            % find possible links for this end
            linkPossible = find( uniqueEnds( jLink) == [linkCosts.idx1]);
            
            [~, minIdx] = min( [linkCosts(linkPossible).cost]);
            
           if linkCosts( linkPossible( minIdx) ).completeLink
               
               imLink = 0 * brokenNetwork;
               
               imLink( [ uniqueEnds( jLink),  linkCosts( linkPossible( minIdx) ).idx2] ) = 1;
               
               imLink = bwconvhull( imLink);
               
               finalNetwork = logical( finalNetwork + imLink);
               
           end

        end
        if isfield(linkCosts, 'completeLink')
            successLinks = sum( [ linkCosts.completeLink]);
            rejectLinks = length( linkCosts) - successLinks;
        else
            successLinks = 0;
            rejectLinks = 0;
        end
        fprintf('Successful Links = %d , Rejected Links = %d\n', successLinks, rejectLinks) ;
        
    end

end
%% }}}
