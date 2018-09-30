function [coords, orients] = searchForHighIntensityLines2( imageGaussian, imMask, estParams)

    % To estimate the lines, we'll take the following approach:
    % 1) Start at the brightest point in the image, and find the optimal orientation of any underlying lines  by radially integrating intensity.
    % 2) Once we have a start point and a direction, we continue searching for the underlying line by looking for intensity peaks until a threshold is met.
    % 3) We take this estimated microtubule, and the pixel region it spans and set it to the median value of the original image. We say that we have removed the intensity due to this microtubule.
    % 4) We repeat steps 1-3 until the microtubule estimates are smaller than a certain length
    % 5) We will fit a polynomial through each initial estimate (quadratic or
    %    cubic). This will be our final estimate for the microtubules.

    % extract parameters
    vis = estParams.visibility;
    fov = estParams.fieldOfView;
    step = estParams.stepSize;
    minMTlength = estParams.minMTlength;
    maxIter = estParams.maxMTsearchIterations;
    maxFails = 2; 

    % parse input image and get useful properties 
    imageGaussian = mat2gray( imgaussfilt( imageGaussian, 0.5) );
    imEst = imageGaussian.*imMask;
    imVals = imageGaussian( imMask);

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
    intThreshMTOC = mu+2*sig;
    intMax = 1.2*intThreshMTOC;

    while intMax > intThreshMTOC && iter < maxIter && fails < maxFails
         
        iter = iter +1;

        % Find Bright Point in Image
        [ intMax, idxMax] = max( imEst(:) );
        [ yMax, xMax] = ind2sub( size(imEst), idxMax);

        % check if there is a microtubule at this point by looking for peaks in all directions
        [coordTemp, success, orientStart] = estimateNextPoint( [ xMax; yMax], 0, step, vis, 540, imEst);
        
        
        
        % if a microtubule is found, we treat this point as a possible MTOC
        if success

            angOrient = [ orientStart, orientStart+pi];
            coordCurr = [];
            addedMT = 0;  

            for jAng = 1:2

                ang = angOrient( jAng);
                % Continue finding next points
                coordMT{jAng} = estimateMicrotubulePoints( [xMax; yMax], ang, imEst, step, vis, fov);
                % Only save if mt length exceeds minimum acceptable estimate length and intensity through the micorotubule is bigger than our threshold.
                cdX = []; cdY = [];
                for jC = 1 : size(coordMT{jAng}, 2)-1
                    cdX = round( [ cdX, linspace( coordMT{jAng}(1,jC), coordMT{jAng}(1,jC+1), 10) ]);
                    cdY = round( [ cdY, linspace( coordMT{jAng}(2,jC), coordMT{jAng}(2,jC+1), 10) ]);
                end
                indMT = sub2ind( size(imEst), cdY, cdX ); 
                if size(coordMT{jAng}, 2)*step > minMTlength 
                    addedMT = 1;
                    coords = { coords{:}, coordMT{jAng} };
                    orients = [orients, ang];
                end
            end

            coordCurr = [ coordMT{1}, coordMT{2} ]; 
            im1 = imEst;
            imMT = logical( 0*imEst);
        try
            cdX = []; cdY = [];
            for jC = 1 : size(coordCurr, 2)-1
                cdX = [ cdX, linspace( coordCurr(1,jC), coordCurr(1,jC+1), 10) ];
                cdY = [ cdY, linspace( coordCurr(2,jC), coordCurr(2,jC+1), 10) ];
            end

            idxMT = sub2ind( size(imEst), round(cdY), round( cdX) );
            imMT( idxMT) = 1;
            imMT = imdilate( imMT, strel('disk', 2) );
            idxMT = find( imMT);
            rV = randValueGenerator( muMin, sig/2, size(idxMT), muMin+sig/2, 0 );
            imEst( idxMT) = rV;
            minV = min( imEst(:) ); maxV = max( imEst(:) ); 
            imEst = mat2gray( imEst) * ( maxV-minV) + minV;
            imEst = imEst .* imMask;
        end
            % de-intensity the mtoc point
            spread=4;
            for jX = 1 : size(imEst, 2); for jY = 1 : size(imEst, 1)
                if norm( [jX-xMax, jY-yMax] ) <= spread
                    imEst(jY, jX) = randValueGenerator( muMin, sig/2, 1, muMin+sig/2, 0);
                end
            end; end
            imEst = imEst.*imMask;
            im2 = imEst;
    plotOn = 0;
        if plotOn
            figure('NumberTitle', 'off', 'Name', ['Iter_' num2str(iter)] )
            subplot(221)
            imagesc(imageGaussian.*imMask); colormap gray
            title('Original Image')
            subplot(222)
            imagesc( im1); colormap gray; hold on; plot( xMax, yMax, 'b*', 'MarkerSize', 14); hold off;
            title( 'Search Image: MTOC' )
            subplot(223)
            imagesc( im1); colormap gray; hold on; plot( coordMT{1}(1, :), coordMT{1}(2,:), 'r-', 'LineWidth', 2); plot( coordMT{2}(1, :), coordMT{2}(2,:), 'r-', 'LineWidth', 2); plot( xMax, yMax, 'b*', 'MarkerSize', 14); hold off
            title( 'Microtubule found' )
            subplot(224)
            imagesc(im2); colormap gray; title('Image for next iteration')
        end

        else
            spread=4;
            for jX = 1 : size(imEst, 2); for jY = 1 : size(imEst, 1)
                if norm( [jX-xMax, jY-yMax] ) <= spread
                    imEst(jY, jX) = randValueGenerator( muMin, sig/2, 1, muMin+sig/2, 0);
                end
            end; end
            imEst = imEst.*imMask;
            fails = fails+1;
            disp(fails)
        end


    end

    % as a final step, we say that the brightest point is definitely on a MT, and that only nucleations near that MT are allowed, so we get rid of any MTs that nucleate farther away than the diameter of the nucleas from the brightest point.
    if length(coords) >= 1
        mtoc_Brightest = round( coords{1}(:,1)); xMax = mtoc_Brightest(1); yMax = mtoc_Brightest(2);
        imMTOC = 0*imageGaussian; imMTOC( sub2ind( size(imageGaussian), yMax, xMax ) ) = 1;
        spread=40;
        for jX = 1 : size(imageGaussian, 2); for jY = 1 : size(imageGaussian, 1)
            if norm( [jX-xMax, jY-yMax] ) <= spread
                imMTOC(jY, jX) = 1;
            end
        end; end
        
        % check which mtocs are in this region and only keep those MTs
        angsRM = [];
        for jAng = 1 : length( coords)
            idxMT = sub2ind( size(imageGaussian), coords{jAng}(2,1), coords{jAng}(1,1) );
            if imMTOC( idxMT) == 0
                angsRM = [angsRM, jAng];
            end
        end
        disp(angsRM)
        coords( angsRM) = [];
        orients( angsRM) = [];
    end



    function rVals = randValueGenerator( mu, sig, sizeRVal, ub, lb)
    % randValueGenerator: generates a random value from a normal distribution centered at mean mu with standard deviation sig
        rVals = mu + ( sig * randn( sizeRVal) ); 
        rVals( rVals > ub) = ub;
        rVals( rVals < lb) = lb;
    
    end

end

