function networkout = destroyNetworkJunctions( networkIn)
        
        % bomb places where their are any connections
        
        % clear the border
        networkIn( 1, :) = 0; networkIn( :, 1) = 0; networkIn( end, :) = 0; networkIn( :, end) = 0;
        networkout = networkIn;
        cc = bwconncomp( networkout); nObj = cc.NumObjects;
        countBombed1 = 0; countBombed2 = 0; countBombed3=0;
        for jObj = 1 : nObj

            % get all its pixels, and turn on just the object
            pixList = cc.PixelIdxList{ jObj};
            imObj = 0 * networkout;
            imObj( pixList) = 1;
            % get cartesian indices of pixel
            [yList, xList] = ind2sub( size(imObj), pixList );

            % check if any pixel satisfies our criterion for being a connection
            % pixel
            for jPix = 1 : length(pixList)
                if pixList( jPix) == 11400
                    stopH = 1;
                end
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
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
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
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
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
                            ( ( sum( Mindist == 1) == 4) && ( sum( Mindist >= 1.3) == 1) ) ||...
                            ( sum( Mindist >= 2) == 2)
                        [y,x] = ind2sub( size(networkout), pixList( jPix) );
                        networkout( y-1:y+1, x-1:x+1 ) = 0;
                        countBombed3 = countBombed3+1;
                    end
                elseif sum( nhood(:))  >= 6
                    [y,x] = ind2sub( size(networkout), pixList( jPix) );
                    networkout( y-1:y+1, x-1:x+1 ) = 0;
                    countBombed3 = countBombed3+1;
                end        

            end

        end
        
end


