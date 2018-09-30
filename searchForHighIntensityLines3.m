function [coords, orients] = searchForHighIntensityLines3( imGauss, imLines, imMask, estParams)

    % To estimate the lines, we'll take the following approach:
    % 1) Start at the brightest point in the image, and find the optimal orientation of any underlying lines  by radially integrating intensity.
    % 2) Once we have a start point and a direction, we continue searching for the underlying line by looking for intensity peaks until a threshold is met.
    % 3) We take this estimated microtubule, and the pixel region it spans and set it to the median value of the original image. We say that we have removed the intensity due to this microtubule.
    % 4) We repeat steps 1-3 until the microtubule estimates are smaller than a certain length
    % 5) We will fit a polynomial through each initial estimate (quadratic or
    %    cubic). This will be our final estimate for the microtubules.

    plotflag_iter = 0;
    flag_closeNucleation = 0;
    mtocSpread = 40; % only used when flag_closeNucleation is active

    % extract parameters
    vis = estParams.visibility;
    fov = estParams.fieldOfView;
    step = estParams.stepSize;
    minMTlength = estParams.minMTlength;
    maxIter = estParams.maxMTsearchIterations;
    maxFails = 2; 
    imEst = imLines.*imMask;

    % Get pixel intensity distribution information from the intensity image
    imVals = imGauss( imMask);
    med = median( imVals );
    otsu = multithresh( imVals, 1);
    mu = mean( imVals);
    sig = std( imVals);
    muMin = min( [mu, med]);

    disp( sprintf( 'median=%.3f', med ) ) 
    disp( sprintf( 'mean=%.3f', mu ) )
    disp( sprintf( 'std=%.3f', sig) )
    disp( sprintf( 'otsu=%.3f', otsu(1) )  )
    
    % initialize parameters for looping
    fails = 0;
    iter = 0;
    coords = {};
    orients = [];
    intThreshMTOC = otsu;
    intMax = 1.2*intThreshMTOC;

    while (intMax > intThreshMTOC || max( logical(imEst(:)).*imGauss(:)) > intThreshMTOC) && iter < maxIter && fails < maxFails
         
        iter = iter +1;
        
        % Find Bright Point in Image
        [ intMax, idxMax] = max( imEst(:) );
        [ yMax, xMax] = ind2sub( size(imEst), idxMax);

        % Does MT pass through this Bright Point? 
        [~, success, orientStart] = estimateNextPoint( [ xMax; yMax], 0, step, vis, 540, imEst);
        
        % if MT passes, search for MT on both sides 
        if success
            
            % 2 angles because 2 sides
            angOrient = [ orientStart, orientStart+pi]; 
            for jAng = 1:2
                ang = angOrient( jAng);
                coordMT{jAng} = estimateMicrotubulePoints( [xMax; yMax], ang, imEst, step, vis, fov);
            end

            % find real MTOC using the intensity image informaiton instead of the steerable filtered image
            [coordMT, orientMT, intMax, mtoc] = findRealMTOCViaIntensityImage( coordMT, imGauss);
            intMax

            % append found MTs if min length bigger than some threshold and intensityMTOC is biggger than threshold
            [coords, orients, addedMT] = appendMTifConditionsMet( coords, coordMT, orients, orientMT, minMTlength, step, intMax, intThreshMTOC);

            % Remove the MT from the image, and de-gauss the steetable MTOC
            imPost = removeMicrotubuleFromImage( imEst, coordMT, [xMax, yMax], 3, 0, 0, imMask);
            
            % display iteration plot if flag activated
            displayIterationPlot( imGauss, imEst, imPost, imMask, coordMT, [xMax, yMax], mtoc, iter, plotflag_iter); 
            
        else
            % de-gauss the steetable MTOC
            imPost = removeMicrotubuleFromImage( imEst, [], [xMax, yMax], 3, 0, 0, imMask);
            fails = fails+1;
            disp( sprintf( 'iter = %d , fails = %d', iter, fails) )
        end
        imEst = imPost;
    end
    
    % as a final step, we say that the brightest point is definitely on a MT, and that only nucleations near that MT are allowed, so we get rid of any MTs that nucleate farther away than the diameter of the nucleas from the brightest point.
    [coords, orients] = keepMTOCScloseTogether( coords, orients, imGauss, imMask, mtocSpread, flag_closeNucleation);

     coords = connectingSimilarMicrotubules( coords, imGauss); 


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
    function [coords, orients, addedMT] = appendMTifConditionsMet( coords, coordMT, orients, orientMT, minMTlength, step, intMax, intThresh)
        
        addedMT = 0;
        for jAng = 1 : 2
            if size( coordMT{jAng}, 2)*step > minMTlength && intMax >= intThresh
                addedMT = 1;
                coords = { coords{:}, coordMT{jAng} };
                orients = [orients, orientMT(jAng)];
            end
        end

    end
    %  }}}
    
    % removeMicrotubuleFromImage {{{
    function imPost = removeMicrotubuleFromImage( imPrior, coordMT, deGaussPoint, deGaussPointSpread, mtVal, ptVal, imMask)
        
        imPost = imPrior;
        if isempty( coordMT), skipMTremoval = 1; else, skipMTremoval = 0; end
        
        if ~skipMTremoval
            coordCurr = [ coordMT{1}, coordMT{2} ]; 
            imMT = logical( 0*imPrior);
            try
                cdX = []; cdY = [];
                for jC = 1 : size(coordCurr, 2)-1
                    cdX = [ cdX, linspace( coordCurr(1,jC), coordCurr(1,jC+1), 10) ];
                    cdY = [ cdY, linspace( coordCurr(2,jC), coordCurr(2,jC+1), 10) ];
                end

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
        imagesc( imPrior); colormap gray; hold on; plot( coordMT{1}(1, :), coordMT{1}(2,:), 'r-', 'LineWidth', 2); plot( coordMT{2}(1, :), coordMT{2}(2,:), 'r-', 'LineWidth', 2); plot( mtocReal(1), mtocReal(2), 'b*', 'MarkerSize', 14); hold off
        title( 'Microtubule found' )
        subplot(224)
        imagesc(imPost); colormap gray; title('Image for next iteration')

    end
    % }}}

    % keepMTOCScloseTogether {{{
    function [coords, orients] = keepMTOCScloseTogether( coords, orients, imOrg, imMask, spread, runflag) 
        
        if ~runflag, return, end

        if length(coords) >= 1

            mtoc_Brightest = round( coords{1}(:,1)); 
            xMax = mtoc_Brightest(1); yMax = mtoc_Brightest(2);
            imMTOC = 0*imOrg; imMTOC( sub2ind( size(imOrg), yMax, xMax ) ) = 1;
            
            for jX = 1 : size(imOrg, 2); for jY = 1 : size(imOrg, 1)
                if norm( [jX-xMax, jY-yMax] ) <= spread
                    imMTOC(jY, jX) = 1;
                end
            end; end
            
            dispImg( imMTOC.*imOrg.*imMask); hold on;
            for jC = 1 : length(coords)
                plot( coords{jC}(1,1), coords{jC}(2,1), 'r*', 'MarkerSize', 12)
            end, hold off

            % check which mtocs are in this region and only keep those MTs
            angsRM = [];
            for jAng = 1 : length( coords)
                idxMT = sub2ind( size(imGauss), coords{jAng}(2,1), coords{jAng}(1,1) );
                if imMTOC( idxMT) == 0
                    angsRM = [angsRM, jAng];
                end
            end
            disp(['MT removed: idx = ', num2str(angsRM) ]);
            coords( angsRM) = [];
            orients( angsRM) = [];
        end

    end
    % }}}

    % connectingSimilarMIcrotubules {{{
    function coords = connectingSimilarMicrotubules( coords, imInt) 

        % Every microtubule can be linked to every other microtubule
        % A Microtubule cannot be linked to its mtoc pair which nucleates on the other side
        
        linkAccept = 15;
        phiScale = 300;
        intScale = 10;
        lenScale = 1;
        numMT = length(coords);

        linkCosts = zeros( numMT);
        linkType = cell( numMT);
        % Analyze each Microtubule/Path and get its properties
        
        for jmt = 1 : numMT
            lineProps( jmt) = analyzeLineImage( imInt, coords{jmt});
        end
        for jmt = 1 : numMT
            [lineProps( jmt).contourLength, lineProps(jmt).endLength] = findContourLength( coords{jmt});
        end

        for j1 =  1: numMT, for j2 = 1: numMT
           
            %  dont connect if same MT
            if j1 == j2; continue; end

            %  dont conenct if other mt has same mtoc as current mt
            if all( coords{j1}(:,1) == coords{j2}(:,1) ), continue; end

            m1 = lineProps(j1); m2 = lineProps(j2);

            % cost depends on three quantities
            %  1) link length
            %  2) link Angle with mt1
            %  3) link ANgle with mt2
            %  4) link intensity vs obj intensities

            % end of mt1 can be linked to either the mtoc of mt2 or to the end of mt2
            end1 = coords{j1}(:, end);
            end2_1 = coords{j2}(:,1);
            end2_2 = coords{j2}(:,end);

            % possible link is to whatever is closest
            if norm( end1 - end2_1) < norm( end1 - end2_2)
                end2 = end2_1; linkObj = 'start';
            else, end2 = end2_2; linkObj = 'tip'; end
            linkType{j1, j2} = linkObj;

            % Generate link and get its properties
            link_length = norm( end1 - end2);
            costLen = lenScale*link_length;

            % find effect of adding link on contour length / end length
            if strcmp( linkObj, 'start')
                [CLnew, ELnew] = findContourLength( [coords{j1}, coords{j2}]);
            elseif strcmp( linkObj, 'tip')
                [CLnew, ELnew] = findContourLength( [coords{j1}, flip(coords{j2}, 2)]);
            end
            energyNew = CLnew/ELnew;
            % check effect of linking on mt 1
            cl1 = m1.contourLength; cl2 = m2.contourLength;
            el1 = m1.endLength; el2 = m2.endLength;
            energyOld = mean( [cl1/el1 cl2/el2]);
            dE = (energyNew - energyOld)/energyOld;
            costE = phiScale*dE;

            % find intensities of mts and compare to intensity of link
            int1 = findUnderlyingIntensity( coords{j1}, imInt);
            int2 = findUnderlyingIntensity( coords{j2}, imInt);
            intLink = findUnderlyingIntensity( [end1 , end2], imInt);
            dInt = abs(intLink - mean([int1, int2]) )/ mean([int1, int2]);
            costInt = intScale*dInt; 

%             disp( sprintf('Possible Link: %d and %d', j1, j2))
%             disp( sprintf('Link Length = %.1f', link_length))
%             disp( sprintf('Cont: old = %.1f, new = %.1f', mean( [cl1, cl2]),  CLnew))
%             disp( sprintf('End: old = %.1f, new = %.1f', mean([ el1, el2 ]), ELnew))
%             disp( sprintf('Energy: old = %.1f, new = %.1f', energyOld, energyNew))
             
            linkCosts(j1, j2) = costLen + costE + costInt;

%             if linkCosts(j1, j2) <linkAccept 
%                 disp('-------------------------')
%                 disp(['Link Successful: ' num2str(j1) ' to ' num2str(j2)])
%                 disp( ['costPhi = ' num2str(costE) ])
%                 disp(['costLen = ' num2str(costLen)])
%                 disp(['costInt = ' num2str(costInt)])
%                 disp(['linkCost = ' num2str(linkCosts(j1, j2)) ])
%                 disp('-------------------------')
%             end

        end, end

        linkCosts( linkCosts == 0) = max( linkCosts(:) );
        
        % return the links in pairs
        indLink = find( linkCosts < linkAccept);
        [mt1, mt2] = ind2sub( size(linkCosts), indLink);
        
        % scrollthrough these and pick the best ones (for pairs pick the one with lowest cost)
        links = {}; costs=[];
        for j1 = 1 : length(mt1); 

            minCost = linkCosts( mt1(j1), mt2(j1) );
            match1 = mt1(j1);
            match2 = mt2(j1);

            for j2 = 1 : length(mt2)
                c1 = linkCosts( mt1(j1), mt2(j1) );
                c2 = linkCosts( mt1(j2), mt2(j2) );
                %if multiple links per 1st mt, only pick best
                if mt1(j1) == mt1(j2) && c2  < minCost
                    minCost = c2;
                    match2 = mt2(j2);
                end
            end

            % check if opposite is opposite and which is best match
            for j2 = 1 :length(mt2)
                c1 = linkCosts( mt1(j1), mt2(j1) );
                c2 = linkCosts( mt1(j2), mt2(j2) );
                if mt1(j1) == mt2(j2) && mt2(j1) == mt1(j2) && c2 < minCost
                    match1 = mt1(j2);
                    match2 = mt2(j2);
                    minCost = c2;
                end
            end

            if any(cellfun( @(c) ((c(1)==match1 && c(2)==match2) || (c(1)==match2 && c(2)==match1)), links));
                continue
            end
            links = { links{:}, [match1, match2]};
            costs = [costs, minCost];
        end
        
        indRM = [];
        for jDis = 1 : length( costs)
            
            m1 = links{jDis}(1); m2 = links{jDis}(2);
            indRM = [indRM, m2]; 
            switch linkType{m1, m2}
                case 'tip'
                coords{m1} = [ coords{m1} , coords{m2}]; 
                case 'start'
                coords{m1} = [ coords{m1}, flip( coords{m2}, 2) ];
                otherwise
                error('Something went horribly wronggg')
            end

            disp('-------------------------')
            disp(['Link Successful: ' num2str(m1) ' to ' num2str( m2 )])
            disp(['linkCost = ' num2str( costs(jDis)) ])
            disp(['linkType = tip - ' linkType{ m1, m2 }])
            disp('-------------------------')
        
        end

        coords( indRM) = [];

        function int = findUnderlyingIntensity( coords, imInt)

            x = coords(1,:); y = coords(2,:);
            int = [];
            for jx = 1:length(x)-1
                im0 = 0*imInt;
                im0( y(jx), x(jx) ) = 1;
                im0( y(jx+1),x(jx+1) ) = 1;
                im0 = bwconvhull( im0);
                int = [int; imInt(find(im0))];
            end
            int = mean(int);

        end

    end
    % }}}

    %  findContourLength {{{
    function [contourLength, endLength] = findContourLength( coords)

        x = coords(1,:);
        y = coords(2,:);
        contourLength = sum( sqrt(diff(x).^2 + diff(y).^2) );
    
        endLength = norm( coords(:,1) - coords(:, end) );

    end
    %  }}}

    % shapeFilter {{{
    function imFilt = shapeFilter( imPlane, minSolidity)

        cc = bwconncomp( logical(imPlane));
        st = regionprops( cc, 'Solidity');
        disp( [st(:).Solidity])

        % keep objects whose solidity is greater than 0.5
        idxKeep = find( [st(:).Solidity] > minSolidity);
        imFilt = 0*imPlane;
        for jKeep = 1 : length(idxKeep)
            imFilt( cc.PixelIdxList{idxKeep(jKeep)} ) = 1;
        end
        imFilt = imFilt .* imPlane;

    end
    % }}}

end
