function isolatedComponents = findMinimalComponents( newNetwork)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

plotFlag = 1;

imBombed = bombNetworkJunctions( newNetwork);
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
    [ fullSeq, runFlag] = extendSequences( endPoints, fullSeq, initSeq);
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

% all plots
if plotFlag
    
    %     dispImg( newNetwork, imBombed, [1 2])

    %     imEnd = 0 * imBombed;
    %     imEnd( [endPoints.pixId] ) = 1;
    %     imEnd2 = zeros( [ size(imEnd), 3] );
    %     imEnd2(:,:,1) = imBombed-imEnd;
    %     imEnd2(:,:,2) = imEnd;
    %     dispImg( imEnd2); title('Endpoints of Bombed Network'); set(gca, 'FontSize', 15);

    nR = 2; nC = ceil( length( fullSeq)/2);
    dispImg( imConnections{:}, [nR, nC] )
    
end






isolatedComponents = 1;


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
        
        % for each object in the broken network
        for jObj = 1 : nObj
            
            % get all its pixels, and turn on just the object
            pixList = cc.PixelIdxList{ jObj};
            imObj = 0 * imBombed;
            imObj( pixList) = 1;
            % get cartesian indices of pixel
            [yList, xList] = ind2sub( size(imObj), pixList );
            
            % check if any pixel satisfies our criterion for being an
            % endpixel.
            for jPix = 1 : length(pixList)
                
                y = yList( jPix); x = xList( jPix);
                
                % get the local 3x3 neighborhood
                nhood = imObj( y-1:y+1, x-1:x+1);
                
                % apply criterion for endpoint
                if sum( nhood( 2:2:end) ) < 2
                    % if there aren't any horizontal/vertical or L-shaped
                    % connections, it is an endpoint if there is only one
                    % turned on pixel in the neighborhood. It is also an
                    % endpoint if there are 2 turned on pixels but they are
                    % very close together, indicating that something ends
                    % at the current pixel.
                    
                    if sum( nhood(:) ) - 1 == 1
                        
                        endPoints = [ endPoints; pixList( jPix), x, y, jObj];
                        
                    elseif sum( nhood(:) ) - 1 == 2

                        % find the pixels that are turned on
                        idx = find( nhood);
                        idx = idx( idx ~= 5); % not the current pixel
                        [yi, xi] = ind2sub( size(nhood) , idx );

                        % find dist between pixels
                        if norm([diff(xi) diff(yi)]) < 2.1
                            
                            endPoints = [ endPoints; pixList( jPix), x, y, jObj];
                            
                        end
                        
                    end
                    
                end
                
            end

        end
        
        for jEnd = 1 : size( endPoints, 1)
            
            endPointsStruct( jEnd).pixId = endPoints( jEnd, 1);
            endPointsStruct( jEnd).x = endPoints( jEnd, 2);
            endPointsStruct( jEnd).y = endPoints( jEnd, 3);
            endPointsStruct( jEnd).objId = endPoints( jEnd, 4);
            
        end
        
        endPoints = endPointsStruct;
        
    end

    function endPoints = labelEndPoints( endPoints)
        
        for j = 1 : length(endPoints)
            endPoints(j).label = j;
        end
        
    end

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
            seq( isnan( seq) ) = [];
            initSeqCell{ jSeq} = seq;
        end
        initSeq = initSeqCell;
        
    end

    function [fullSeq, extendFlag] = extendSequences( endPoints, initSeq, initPairs)
        % extendSequences : extends sequences by linking endpoints if they
        % are close in space and their orientation matches
        
        phiSpreadCone = deg2rad( 35);
        lengthCone = 10;
        maxPhiDiff = deg2rad(80);
        fullSeq = initSeq;
        extendFlag = 0;
        
        % for each endpoint, find its ideal links
        for jNode = 1 : length( endPoints)
            
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
                imCone( round( node1y + lengthCone * sin( jPhi)), round( node1x + lengthCone * cos( jPhi) ) ) = 1;
            end
            imCone = bwconvhull( imCone);
            
            % Now, generate image containing any objects that appear in the
            % search cone. First Create image with all objects except the
            % one belonging to our endpoint.  
            imOther = imBombed;
            imOther( cc.PixelIdxList{ endPoints( jNode).objId} ) = 0;
            imSearch = imCone .* imOther;
            
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

                    % find difference of phi values to see if there is a match
                    phiDiff = min( abs( (phi2 + phi1) + pi*[-3:2:3] ) );

                    % we need to determine what max linking angle to allow. We
                    % will base this on how far the link point is. 
                    % phiMax( d=1) = 30 degrees, phiMax( d=dMax) = 60 degrees
                    [ linkY, linkX] = ind2sub( size(newNetwork), node2Idx);
                    dist = norm( [ linkY-node1y, linkX-node1x ] );
                    phiMax = maxPhiDiff * ( 1 + 1/(lengthCone-1) * dist);

                    if phiDiff < phiMax
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
                        if cSeq(2) ~= jLink
                            
                            % find link pair
                            pairLink = unique( cellfun( @(v) (v(1)==jLink)*v(2),  initPairs ));
                            pairLink = pairLink( pairLink ~= 0);
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
                            pairLink = pairLink( pairLink ~= 0);
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

    function lineProps = findLineImageVariance( lineImage)
        % findCurveVariance: find the variance of the local curvatures of
        % an object as well as some other properties
        
        % get the ordering of the connections in the line image
%         [B,~] = bwboundaries( lineImage, 4, 'noholes');
%         boundaryPixels = B{1};
%         numPix = length( boundaryPixels);

        % find endpoint of line
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
        filtSize = 5;
        
        thetaList = [];
        for jPix = 1 : numPix
            
            if jPix <= (filtSize-1)/2
                startPix = 1;
            else
                startPix = jPix - (filtSize-1)/2;
            end
            if jPix > numPix - (filtSize-1)/2
                endPix = numPix;
            else
                endPix = jPix + (filtSize-1)/2;
            end
            
            localLine = xyList( startPix : endPix, :)';
            
            xD = diff(localLine(1,:) );
            yD = diff(localLine(2,:) );
            
            thetaList = [ thetaList, mean( atan2(yD, xD) ) ];
            
        end
        
        % store information about line/curve
        lineProps.length = numPix;
        lineProps.theta = mean( thetaList);
        lineProps.thetaVar = std( thetaList);
        lineProps.thetaRaw = thetaList;
        
    end

    function lineProps = analyzeLineImages( imageCellArray)
        % analyzeLineImages: runs the findLineImageVariance function on
        % each image inside the imageCellArray and stores in the lineProps
        % structure
        
        for jImg = 1 : length( imageCellArray)
            
            lineProps(jImg) = findLineImageVariance( imageCellArray{ jImg} );
            
        end
        
    end

    function BestRoads= findBestGlobalSequences( allRoads, lineProps)
        % findBestGlobalSequences: finds the minimal sequences that best
        % explains the original image by enforcing decisions on strongly
        % overlapping roads.
        
        maxOverlap = 15;
        
        numRoads = length( allRoads);
        
        winners = []; losers = [];
        for jRoad = 1 : numRoads
            
            % get overlaps
            overlaps = cellfun( @(v) sum( allRoads{jRoad}(:) & v(:) ) > maxOverlap, allRoads );
            
            % remove overlap with itself
            overlaps(jRoad) = 0;
            
            conflictRoads = find( overlaps);
        
            % for each conflicting road, check the oreintations of all the
            % other segments to see which roadway this overlap truly
            % belongs to.
            winner = jRoad;
            for kRoad = conflictRoads
                
                var1 = lineProps( jRoad).thetaVar;
                var2 = lineProps( kRoad).thetaVar;
                
                % find propereties of overlap region
                imOverlap = logical( 0 * allRoads{1} );
                
                imOverlap( allRoads{jRoad}(:) & allRoads{kRoad}(:) ) = 1;
                overlapProps = findLineImageVariance( imOverlap);
                
                % get local orientations of overlap region
                thBoth = overlapProps.thetaRaw;
                th1 = lineProps(jRoad).thetaRaw;
                th2 = lineProps(kRoad).thetaRaw;
%                 loc1 = strfind( th1, thBoth);
%                 loc2 = strfind( th2, thBoth);
%                 th1( loc1: loc1+length(thBoth)-1) = [];
%                 th2( loc2: loc2+length(thBoth)-1) = [];
                
                th1 = smooth(th1, 10); th2 = smooth(th2, 10); thBoth = smooth( thBoth, 10);
%                 figure; plot(th1); hold on; plot(th2); plot(thBoth);
                
                d1 = diff( th1); d2 = diff( th2); dB = diff( thBoth);
%                 figure; plot(d1); hold on; plot(d2); plot(dB);
                
                N = min( [length(d1), length(d2)] );
                % get local orientations of the lines in competition after
                % excluding the overlap region
                
                [h1, p1, ~] = kstest2( dB, d1 );
                [h2, p2, ~] = kstest2( dB, d2 );

                if p1 >= p2
                    winner = jRoad;
                    losers = [ losers, kRoad];
                elseif p1 < p2
                    winner = kRoad;
                    losers = [ losers, jRoad];
                end
                
            end
            
            winners = [ winners, winner];
        
        end
        
        BestRoads = { allRoads{ unique( winners)} };
        
    end

end

