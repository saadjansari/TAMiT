function BestConnections = findMinimalComponents( newNetwork, params, plotFlag)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

imBombed = bombNetworkJunctions( newNetwork);
imBombed = bwareafilt( imBombed, [ 4 Inf]);
cc = bwconncomp( imBombed); nObj = cc.NumObjects;

% We will label all the objects and their endpoints uniquely.
% find the endpoints
endPoints = findNetworkEndpoints( imBombed);

% Each endpoint is labeled uniquely, from 1 to N, but we assign labels anyway
endPoints = labelEndPoints( endPoints);

% We find the local orientation of each endpoint and store it
endPoints = findEndpointOrientation( endPoints);

% To start, we link endpoints that belong to the same object. If an object
% has more than two endpoints, we use combinations of 2 endpoints (one for
% entry and one for exit from the obejct). If an obejct just has one
% endpoint, we say that it has a secondary invisible endpoint (NaN) that
% cannot be used for entry or exit.
initSeq = linkObjectEndpoints( endPoints);

% Extend sequences until no longer possible
runFlag = 1; fullSeq = initSeq; count = 0;
while runFlag && count < 10
    [ fullSeq, runFlag] = extendSequences( endPoints, fullSeq, initSeq, params.costs);
    fullSeq = removeChildrenSequences( fullSeq);
    count = count+1;
end

% Given complete sequences, generate images of sequences. The image of a
% sequence is created by either using a pre-made object that is a part of
% the original image, or if its a connection via a jump ( it is a straight
% line connecting the two nodes.
imConnections = generateImagesFromSequences( imBombed, fullSeq, initSeq, endPoints);

% Difference sequences can share roadways upto a certain length (say 10
% pixels). If there is a roadway overlap bigger than 10, then we will find
% which road the overlap region truly belongs to. We will dispose of the
% other road.
lineProps = analyzeLineImages( imConnections);

% For each road overlap/conflict, get images without the overlap, and
% analyze the orientations and lengths of the other segments to find the
% best match.
BestConnections = findBestGlobalSequences( imConnections, lineProps);

BestConnections = removeSimilarLines( BestConnections);

% Finally keep any parts of imBombed that are not a part of any of the
% best connections as they are
BestImage = cat(3, BestConnections{:}); BestImage = logical(sum(BestImage, 3) );
cc2 = bwconncomp( heaviside(imBombed-BestImage-0.1) );
for jPiece = 1 : cc2.NumObjects
    
    img = 0*imBombed;
    img( cc2.PixelIdxList{jPiece} ) = 1;
    BestConnections = { BestConnections{:}, img};
    
end

% all plots
if plotFlag
    
%     dispImg( newNetwork, imBombed, [1 2])

%     imEnd = 0 * imBombed;
%     imEnd( [endPoints.pixId] ) = 1;
%     imEnd2 = zeros( [ size(imEnd), 3] );
%     imEnd2(:,:,1) = imBombed-imEnd;
%     imEnd2(:,:,2) = imEnd;
%     dispImg( imEnd2); title('Endpoints of Bombed Network'); set(gca, 'FontSize', 15);
% 
%     nR = 2; nC = ceil( length( fullSeq)/2);
%     dispImg( imConnections{:}, [nR, nC] ); title('All Possible Components');
%     
%     nR = 2; nC = ceil( length( BestConnections)/2);
%     dispImg( BestConnections{:}, [nR, nC] ); title('Best Components');
   
    plotFinalComponents( newNetwork, BestConnections);
    
end

%% imBombed = bombNetworkJunctions( newNetwork) {{{
    function imBombed = bombNetworkJunctions( newNetwork)
        
        
        % bomb places where their are any connections
        imBombed = newNetwork;
        cc = bwconncomp( imBombed); nObj = cc.NumObjects;
        countBombed1 = 0; countBombed2 = 0; countBombed3=0;
        for jObj = 1 : nObj

            % get all its pixels, and turn on just the object
            pixList = cc.PixelIdxList{ jObj};
            imObj = 0 * imBombed;
            imObj( pixList) = 1;
            % get cartesian indices of pixel
            [yList, xList] = ind2sub( size(imObj), pixList );

            % check if any pixel satisfies our criterion for being a connection
            % pixel
            for jPix = 1 : length(pixList)

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
                        [y,x] = ind2sub( size(imBombed), pixList( jPix) );
                        imBombed( y-1:y+1, x-1:x+1 ) = 0;
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
                        [y,x] = ind2sub( size(imBombed), pixList( jPix) );
                        imBombed( y-1:y+1, x-1:x+1 ) = 0;
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
                            ( ( sum( Mindist == 1) == 4) && ( sum( Mindist >= 2) == 1) ) ||...
                            ( sum( Mindist >= 2) == 2)
                        [y,x] = ind2sub( size(imBombed), pixList( jPix) );
                        imBombed( y-1:y+1, x-1:x+1 ) = 0;
                        countBombed3 = countBombed3+1;
                    end
                elseif sum( nhood(:))  >= 6
                    [y,x] = ind2sub( size(imBombed), pixList( jPix) );
                    imBombed( y-1:y+1, x-1:x+1 ) = 0;
                    countBombed3 = countBombed3+1;
                end        

            end

        end
        
    end
%% }}}

%%  endPoints = findNetworkEndpoints( imBombed) {{{
    function endPoints = findNetworkEndpoints( imBombed)
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
%             imObj = 0 * imBombed;
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

        imEndPoints = bwmorph( imBombed, 'endpoints');
        
        % for every connected pixels, find its centroid and that will be
        % the endpoint
        ccc = regionprops( imEndPoints, 'Centroid', 'PixelIdxList');
        endP = cat( 1, ccc(:).PixelIdxList);
        for jEnd = 1 : length(endP)
            
            endPointsStruct( jEnd).pixId = endP(jEnd);
            [ endPointsStruct( jEnd).y, endPointsStruct( jEnd).x] = ind2sub( size(imEndPoints), endP(jEnd));
            
            % get the object number associated with this endpoint
            objs = cellfun( @(v) any( v(:) == endPointsStruct( jEnd).pixId), cc.PixelIdxList );
            endPointsStruct( jEnd).objId = find( objs);
            if isempty( find( objs) )
                error('Oh why o why wont...you fix me')
            end
        end
        
        endPoints = endPointsStruct;
       
        % ensure that a single object has a maximum of 2 endpoints. bwmorph
        % has its drawbacks. if more than 2 endpoints, pick the two points
        % that are the farthest.
        objId = [endPoints.objId];
        numObj = max(objId);
        idxKeep = [];
        for jObj = 1 : numObj
            
            endIdx = find(objId == jObj);
            x = [endPoints( endIdx).x];
            y = [endPoints( endIdx).y];
            
            distCurr = 0;
            idx1Keep = endIdx(1);
            idx2Keep = endIdx(2);
            for p1 = 1 : length(x)
                for p2 = 1 : length(x)
            
                    dist = norm( [x(p1)-x(p2), y(p1)-y(p2)] );
                    if dist > distCurr
                        idx1Keep = p1;
                        idx2Keep = p2;
                        distCurr = dist;
                    end
                    
                end
            end
            idxKeep = [ idxKeep, endIdx(idx1Keep), endIdx(idx2Keep)];
        end
        idxRm = setdiff( 1:length(ccc), idxKeep);
        endPoints(idxRm) = [];
    end
%% }}}

%% endPoints = labelEndPoints( endPoints) {{{
    function endPoints = labelEndPoints( endPoints)
        
        for j = 1 : length(endPoints)
            endPoints(j).label = j;
        end
        
    end
%% }}}

%% endPoints = findEndpointOrientation( endPoints) {{{
    function endPoints = findEndpointOrientation( endPoints)
        % findEndpointOrientation : finds the local orientation of
        % endPoints by looking at the closest pixels connected to the
        % endpoint.
        
        % create a solid circle of radius 5 pixels that will be centered at
        % each endpoint to find the closest pixels.
        searchRad = 15;
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
            imObj = 0 * newNetwork;
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

%% initSeq = linkObjectEndpoints( endPoints) {{{
    function initSeq = linkObjectEndpoints( endPoints)
        % linkObjectEndpoints: links endpoints that contain the same object
        % label in all double permutations
        
        % get object indices and the total number of objects
        objIdx = [ endPoints.objId ];
        objIdxUnique = unique( objIdx);
        numObjects = length( objIdxUnique);
        
        % for each object, search for endpoints associated with it, and
        % create combinations of two of the elements
        initSeq = [];
        for jObj = 1 : numObjects
            
            cObj = objIdxUnique( jObj);
            
            % get all labels associated with this object
            objLabels = [endPoints( find( objIdx == cObj) ).label ];
            
            % get all possible forward and reverse combinations of these
            % labels chosen two at a time
            if length( objLabels) >= 2
                initSeq = [ initSeq ; nchoosek( objLabels, 2) ; flip(nchoosek( objLabels, 2), 2) ];
            else 
                initSeq = [ initSeq ; objLabels NaN];
            end
        end
        % all possible 2-choose sequences are now in the rows of the
        % numeric array comb. We will create cell array to store all the
        % sequences
        initSeqCell = cell(1, size(initSeq, 1) );
        for jSeq = 1 : size(initSeq, 1)
            
            seq = initSeq( jSeq, :);
%             seq( isnan( seq) ) = [];
            initSeqCell{ jSeq} = seq;
        end
        initSeq = initSeqCell;
        
    end
%% }}}

%% [fullSeq, extendFlag] = extendSequences( endPoints, initSeq, initPairs, cost) {{{
    function [fullSeq, extendFlag] = extendSequences( endPoints, initSeq, initPairs, cost)
        % extendSequences : extends sequences by linking endpoints if they
        % are close in space and their orientation matches
        
%         phiSpreadCone = deg2rad( 35);
%         lengthCone = 10;
%         maxPhiDiff = deg2rad(80);
        phiSpreadCone = cost.maxPhiDiff_EE/2;
        lengthCone = cost.maxDistLink;
        maxPhiDiff = deg2rad(20); % this is the minimum bound for difference in phis. max is 3 times this at maximum length
        fullSeq = initSeq;
        extendFlag = 0;
        
        % for each endpoint, find its ideal links
        for jNode = 1 : length( endPoints)
            
            if jNode == 16
                stopHe = 1;
            end
            % get node label, node orientation, node index, x and y
            node1 = endPoints( jNode).label;
            phi1 = endPoints( jNode).phi;
            node1Idx = endPoints( jNode).pixId;
            node1x = endPoints( jNode).x;
            node1y = endPoints( jNode).y;
            
            % draw a cone centered at the current node and search for other
            % nodes inside cone
            % Create a mask with just the search cone region turned on
            phiSweep = linspace( phi1-phiSpreadCone, phi1+phiSpreadCone, 50);
            imCone = 0 * newNetwork; 
            imCone( node1y, node1x) = 1;
            for jPhi = phiSweep
                xV = round( round( node1x + lengthCone * cos( jPhi) ) ); xV(xV < 1) = 1; xV( xV > size(imCone,2) ) = size(imCone,2);
                yV = round( node1y + lengthCone * sin( jPhi)); yV(yV < 1) = 1; yV( yV > size(imCone,1) ) = size(imCone,1);
                imCone( yV, xV ) = 1;
            end
            imCone = bwconvhull( imCone);
            
            % Now, generate image containing any objects that appear in the
            % search cone. First Create image with all objects except the
            % one belonging to our endpoint.  
            imOther = imBombed;
            imOther( cc.PixelIdxList{ endPoints( jNode).objId} ) = 0;
            try;imSearch = imCone .* imOther;catch
                err=1;
            end
            
            % We will get the pixel Indices of objects in the cone and
            % compare it with pixel indices of the newNetwork image to
            % determine which objects these belong to
            ccLinks = bwconncomp( imSearch);
            numLinks = ccLinks.NumObjects;
            possibleLinks = [];
            for jLink = 1 : numLinks
                
                % get pixel ids of this object
                pixList = ccLinks.PixelIdxList{ jLink}';
                
                % compare with endpoint indices to determine which endpoint
                % is associated with this object
                nodeLink = find( ismember( [endPoints.pixId], pixList ) );
                
                if ~isempty(nodeLink)
                    
                    if length( nodeLink) > 1
                        % pick the closest one
                        xIdx = [ endPoints( nodeLink).x];
                        yIdx = [ endPoints( nodeLink).y];
                        dist = 0*nodeLink;
                        for jj = 1 : length(nodeLink)
                            dist(jj) = norm( [ yIdx(jj) - node1y , xIdx( jj) - node1x]);
                        end
                        [~, jIdx] = min(dist);
                        nodeLink = nodeLink( jIdx);
%                         error( 'multiple endpoints from the same object detected in the search region');
                    end
                    node2Idx = endPoints( nodeLink).pixId;
                    phi2 = endPoints( nodeLink).phi;
                    node2 = endPoints( nodeLink).label;
                    
                    % determine cost of the link by taking two vectors of
                    % unit unit at phi1 and phi2. In a perfect situation,
                    % the length of v1+v2 will be 0. If phi1 and phi2 are
                    % off by 60 degrees, the length of v1+v2 will be equal
                    % to 1 (equilateral triangle)
                    v1 = [ cos(phi1), sin(phi1) ];
                    v2 = [ cos(phi2), sin(phi2) ];
                    
                    % we need to determine what max linking angle to allow. We
                    % will base this on how far the link point is. 
                    % phiMax( d=1) = 30 degrees, phiMax( d=dMax) = 60 degrees
                    [ linkY, linkX] = ind2sub( size(newNetwork), node2Idx);
                    dist = norm( [ linkY-node1y, linkX-node1x ] );
                    phiMax = maxPhiDiff * ( 1 + 2/(lengthCone) * dist);
                    vMaxAllowed = [ cos(0)+cos( pi/2 + phiMax), sin(0)+sin( pi/2 + phiMax ) ];
                    vMaxAllowed = norm( vMaxAllowed);

                    if norm( v1+v2) < vMaxAllowed
                        possibleLinks = [possibleLinks, node2];
                    end
                    
                end
            end
            
            % Now take any sequences that end with the label for the
            % current node, and append an extra label from possible links.
            % Also take sequences that start with the label for the current
            % node and prepend an extra label from the possible links
            if ~isempty( possibleLinks)
                
                % get first and last entries
                firstEntries = cellfun( @(v) v(1), fullSeq);
                
                % find indices in first and last entries that match the
                % current node label
                startMatch = find(firstEntries == node1);
                rmList = [];
                for jMatch = startMatch
                    
                    % obtain the sequence and remove it from the sequence
                    % list
                    cSeq = fullSeq{ jMatch};
                    
                    % prepend possible links, and store it back in the
                    % sequence list
                    prependFlag = 0;
                    
                    for jLink = possibleLinks
                        if length( cSeq) == 1 || cSeq(2) ~= jLink
                            
                            % find link pair
                            pairLink = unique( cellfun( @(v) (v(1)==jLink)*v(2),  initPairs ));
                            pairLink = pairLink( pairLink ~= 0 & ~isnan( pairLink) );
                            if length(pairLink) > 1
                                error( 'issue with pairLinklength')
                            elseif isempty( pairLink)
                                pairLink = NaN;
                            end
                            fullSeq = { fullSeq{:} , [ pairLink, jLink, cSeq] };
                            extendFlag = 1;
                            prependFlag = 1;
                        end
                    end
                    if prependFlag; rmList = [ rmList, jMatch]; end
                end
                try
                fullSeq( rmList) = [];
                catch
                    err = 1;
                end
                % repeat for end matches.
                lastEntries = cellfun( @(x) x(end), fullSeq);
                endMatch = find(lastEntries == node1);
                rmList = [];
                
                for jMatch = endMatch
                    
                    % obtain the sequence and remove it from the sequence
                    % list
                    try
                    cSeq = fullSeq{ jMatch};
                    catch
                        err = 1;
                    end
                    % append possible links, and store it back in the
                    % sequence list
                    appendFlag = 0;
                    
                    for jLink = possibleLinks
                        if cSeq(end-1) ~= jLink
                            
                            % find link pair
                            pairLink = unique( cellfun( @(v) (v(1)==jLink)*v(2),  initPairs ));
                            pairLink = pairLink( pairLink ~= 0 & ~isnan( pairLink) );
                            if length(pairLink) > 1
                                error( 'issue with pairLinklength')
                            elseif isempty( pairLink)
                                pairLink = NaN;
                            end
                            
                            fullSeq = { fullSeq{:} , [ cSeq, jLink, pairLink] };
                            extendFlag = 1;
                            appendFlag = 1;
                        end
                    end
                    if appendFlag; rmList = [ rmList, jMatch]; end
                end
                try
                    fullSeq( rmList) = [];
                catch
                    err = 1;
                end
            end
        end
        
    end
%% }}}

%% parentSeq = removeChildrenSequences( allSeq) {{{
    function parentSeq = removeChildrenSequences( allSeq)
        % removeChildrenSequences: looks at sequences and determines if a
        % partial sequence (children) appears in a larger sequence
        % (parents). It removes all children sequences
        
        rmList = [];
        for jSeq = 1 : length( allSeq)
            
            % obtain possible child sequence
            cSeq = allSeq{jSeq};
            
            % check all other sequences for existence of parents (as is and
            % as flipped)
            otherSeq = allSeq( [1: jSeq-1, jSeq+1:end]);
            parentCheckFwd = cellfun( @(v) any( strfind( v, cSeq) & length(v)>length(cSeq) ) , otherSeq);
            parentCheckBwd = cellfun( @(v) any( strfind( v, flip(cSeq) ) & length(v)>length(cSeq) ) , otherSeq);
            
            if sum( parentCheckFwd(:) )+sum( parentCheckBwd(:) ) > 0
                rmList = [ rmList, jSeq];
            end

        end
        
        % remove children sequences
        parentSeq = allSeq;
        parentSeq( rmList) = [];
        rmList = [];
        % Also remove any duplicated reverse sequences
        skipList = [];
        for jSeq = 1 : length( parentSeq)
            
            if ~any( jSeq == skipList)
                
                % obtain possible child sequence
                cSeq = parentSeq{jSeq};

                % check all other sequences for existence of parents (as is and
                % as flipped)
                doubleCheck = cellfun( @(v) any( strfind( v, flip(cSeq)) ) && length(v)==length(cSeq) , parentSeq);
                doubleIdx = find( doubleCheck);
                skipList = [ skipList, doubleIdx];
                if ~isempty(doubleIdx)
                    rmList = [ rmList, jSeq, doubleIdx(2:end)];
                end
            end
        end
        
        % remove children sequences
        parentSeq( unique(rmList) ) = [];
        
        
    end
%% }}}

%% imConnections = generateImagesFromSequences( imOriginal, fullSeq, initPairs, endPoints) {{{
    function imConnections = generateImagesFromSequences( imOriginal, fullSeq, initPairs, endPoints)
        % generateImagesFromSequences: generates images from sequences.
        % If a link belongs in initPairs, that means its a connection
        % through a pre-existing link in the original image. Other a linear
        % link will be established between the nodes in the sequence
        
        
        cc = bwconncomp( imOriginal);
        
        numSeq = length( fullSeq);
        imConnections = cell( 1, numSeq);
        
        % for each sequence, we'll get consectuve nodes and link them
        for jSeq = 1 : numSeq
            
            imSeq = 0*imOriginal;
            
            % get the sequence and the number of nodes
            cSeq = fullSeq{ jSeq};
            numNodes = length( cSeq);
            
            % link consecutive nodes
            for j1 = 1 : numNodes -1
                
                % get labels for nodes
                n1 = cSeq( j1);
                n2 = cSeq( j1 + 1);
                
                % check if pair belongs to initPairs
                preLinked = cellfun(@(v) any( strfind( v, [n1 n2]) ), initPairs );
                preLinked = any( preLinked);
                
                % if preLinked, then we use the object associated with
                % the nodes
                if preLinked
                    
                    % get the index of this label, and get the
                    % corresponding object index
                    nodeIdx = find( [endPoints.label] == n1);
                    objIdx = endPoints( nodeIdx).objId;
                    
                    % get the object pixels
                    imSeq( cc.PixelIdxList{objIdx} ) = 1;
                    
                else
                    
                    % get the locations of the two endpoints and connect
                    % them with a straight line
                    node1Idx = find( [endPoints.label] == n1);
                    node2Idx = find( [endPoints.label] == n2);
                    loc = [ endPoints( [ node1Idx, node2Idx] ).pixId ];
                    
                    imTemp = 0*imSeq;
                    imTemp( loc ) = 1;
                    imTemp = bwconvhull( imTemp);
                    imSeq = logical( imSeq+imTemp);
                    
                end
                
            
            end
            
            imConnections{ jSeq} = imSeq;

        end
                
    end
%% }}}

%% lineProps = findLineImageVariance( lineImage) {{{
    function lineProps = findLineImageVariance( lineImage)
        % findCurveVariance: find the variance of the local curvatures of
        % an object as well as some other properties
        
        % get the ordering of the connections in the line image
%         [B,~] = bwboundaries( lineImage, 4, 'noholes');
%         boundaryPixels = B{1};
%         numPix = length( boundaryPixels);

        % find endpoint of line
        lineImage = logical( lineImage);
        bw = bwmorph( lineImage, 'endpoints');
        endPtIdx = find(bw == 1);
        [ y1, x1] = ind2sub( size(lineImage), endPtIdx);
        
        % pick an endpoint and measure geodesic distance to each other
        % pixel to get order list of connected pixels
        mat_dist = bwdistgeodesic( lineImage, x1(1), y1(1), 'quasi-euclidean'); %'quasi-euclidean' for 8-connectivity
        comp = find( lineImage);
        comp(:,2) = mat_dist(comp(:,1));
        ordered_list_ind = sortrows(comp,2);
        [ordered_list_sub(:,1) ordered_list_sub(:,2)] = ind2sub( size(lineImage), ordered_list_ind(:,1) ); 
        
        xyList = ordered_list_sub;
        numPix = size( xyList,1);
        % we'll use a moving filter of 5 pixels to find the local
        % orientations
        filtSize = 10;
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
            
%             xD = diff(localLine(1,:) );
%             yD = diff(localLine(2,:) );
            
            xD = localLine(1, end) - localLine(1, 1);
            yD = localLine(2, end) - localLine(2, 1);
            
            thetaNextPossible = mean( atan2(yD, xD) ) + [-2*pi, 0, +2*pi];
            if ~isempty( thetaList)
                [~, thetaNextIdx] = min( abs( thetaNextPossible - thetaList(end) ) );
            else
                thetaNextIdx = 2;
            end
            thetaList = [ thetaList, thetaNextPossible( thetaNextIdx) ];
            
%             thetaList = [ thetaList, mean( atan2(yD, xD) ) ];
            
        end
        
        filtSize = 20;
        filtHalf = (filtSize-1)/2 * mod(filtSize, 2) + filtSize/2 * ~ mod(filtSize, 2);
        contourLengthList = [];
        endLengthList = [];
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
            
            s1 = smooth( thetaList(startPix:endPix) , 7); 

            conLen = sum( sqrt( diff( s1).^2 + 1.^2 ) );
            endLen = sqrt( diff( s1( [1, end]) ).^2 + (length(s1)-1)^2 );
            
%             p1 = xyList( startPix, :);
%             p2 = xyList( endPix, :);
% 
%             conLen = abs( mat_dist( p1(1), p1(2) ) - mat_dist( p2(1), p2(2) ) );
%             endLen = norm( [ p1(1)-p2(1), p1(2)-p2(2) ] );
            
            contourLengthList = [ contourLengthList, conLen];
            endLengthList = [ endLengthList, endLen];
                        
        end
        
        contourLength = max( mat_dist(:) );
        endLength = norm( [ xyList(1,1)-xyList(end,1), xyList(1,2)-xyList(end,2) ] );
        
        % store information about line/curve
        lineProps.numPix = numPix;
        lineProps.theta = mean( thetaList);
        lineProps.thetaVar = std( thetaList);
        lineProps.thetaRaw = thetaList;
        lineProps.meanContPerL = mean( contourLengthList./endLengthList);
        lineProps.contPerL = contourLength/endLength;
        lineProps.length = endLength;
    end
%% }}}

%% lineProps = analyzeLineImages( imageCellArray) {{{
    function lineProps = analyzeLineImages( imageCellArray)
        % analyzeLineImages: runs the findLineImageVariance function on
        % each image inside the imageCellArray and stores in the lineProps
        % structure
        
        for jImg = 1 : length( imageCellArray)
            
            lineProps(jImg) = findLineImageVariance( imageCellArray{ jImg} );
            
        end
        
    end
%% }}}

%% BestRoads= findBestGlobalSequences( allRoads, lineProps) {{{ 
    function BestRoads= findBestGlobalSequences( allRoads, lineProps)
        % findBestGlobalSequences: finds the minimal sequences that best
        % explains the original image by enforcing decisions on strongly
        % overlapping roads.
        
        maxOverlap = 10;
        
        numRoads = length( allRoads);
        
        winners = []; losers = [];
        for jRoad = 1 : numRoads
            
            % get overlaps
            overlaps = cellfun( @(v) sum( allRoads{jRoad}(:) & v(:) )>maxOverlap, allRoads);
            
            % remove overlap with itself
            overlaps(jRoad) = 0;
            
            % find conflict Roads
            conflictRoads = find( overlaps);
            
            if isempty( conflictRoads)
                winners = [winners, jRoad];
            end
            % This road has conflicts with other roads, but there may be
            % different regions of overlap. So we need to get the pixel
            % indices of the overlap regions with each other road and try
            % to group the competition.
            groups = {};
            for iRoad = conflictRoads
                group = [iRoad];
                for iiRoad = conflictRoads
                    if iiRoad ~= iRoad
                        % get pixel indices of both roads
                        r1 = allRoads{ iRoad};
                        r2 = allRoads{ iiRoad};
                        
                        overlapFlag = sum( r1(:) & r2(:))>maxOverlap;
                        if overlapFlag
                            group = [ group, iiRoad];
                        end
                    end
                end
                groups = { groups{:}, group };
            end
            
            % There will be a single winner from each group. At the end,
            % any roads that are not winners will be deemed a loser.
            for jComp = 1 : length( groups)
                cGroup = [ jRoad, groups{ jComp}];
                
                [minVal, winnerIdx] = min([lineProps(cGroup).meanContPerL]);
                
                % if cont lenths too close, pick longer road
                if any([lineProps(cGroup).meanContPerL]-minVal < 1e-4 & [lineProps(cGroup).meanContPerL]-minVal > 0)
                    
                    compete = find( [lineProps(cGroup).meanContPerL]-minVal < 1e-4 & [lineProps(cGroup).meanContPerL]-minVal > 0);
                    competeLabels = cGroup( [ winnerIdx, compete]);
                    [~, winnerIdx] = max( [lineProps( competeLabels ).length] );
                    winnerIdx = find( cGroup == competeLabels(winnerIdx) );
                    
                end
                losers = unique([ losers, setdiff( cGroup, cGroup(winnerIdx) )]);
                winners = unique([ winners, cGroup(winnerIdx)]);
            end
        end
        winnerZ = setdiff( 1:numRoads, losers);
        BestRoads = { allRoads{ winnerZ} };
    end
%% }}}

%% plotFinalComponents( imAll, bestComp) {{{
    function plotFinalComponents( imAll, bestComp)
        
        colors = distinguishable_colors( length(bestComp), {'w', 'k'});
        imBest = 0*imAll; imBest = repmat( imBest, 1, 1, 3);
        
        for jC = 1 : length(bestComp)
            
            [y,x] = find( bestComp{jC});
            
            for jInd = 1 : length(y)
                if any( imBest( y(jInd), x(jInd), :) )
                    imBest( y(jInd), x(jInd), :) = mat2gray( squeeze(imBest( y(jInd), x(jInd), :)) + colors(jC, :)');
                else
                    imBest( y(jInd), x(jInd), :) = colors(jC, :);
                end
            end
            
        end
        dispImg( imAll, imBest, [1 2]);  set( gcf, 'Name', 'network_components', 'NumberTitle', 'off');
%         title( sprintf("Network's minimal components = %d", length(bestComp) ) );
        
    end
%% }}}

%% BestLines = removeSimilarLines( allLines) {{{
    function BestLines = removeSimilarLines( allLines)
        
        % compare each pair of lines to see if any 2 lines are super close
        % if any two lines are super close, keep there average.
        % Define what it means for 2 lines to be similar or super close?
        % if the overlap of any two lines is > 90% of the lines, they are
        % similar
        
        BestLines = {};
        nLines = length( allLines);
        rmVal = [];
        for jLine = 1 : nLines
            if all(jLine ~= rmVal)
                overlaps = cellfun( @(v) sum(allLines{jLine}(:) & v(:)), allLines ); 
                lengs = cellfun(  @(v) sum( v(:) ), allLines );
                ratios = overlaps ./ lengs;

                cLines = find( ratios > 0.9);
                rmVal = unique( [rmVal, cLines( cLines ~= jLine) ] );
                % pick the avg line in cLines
                X = cat(3, allLines{ cLines});
                BestLines = { BestLines{:}, logical( sum( X, 3) ) };
            end
        end
        
    end

end
%% }}}
