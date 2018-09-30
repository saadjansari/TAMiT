function [coords, orients] = searchForHighIntensityLines( imageGaussian, imMask, estParams)

        % To estimate the lines, we'll take the following approach:
        % 1) Start at the brightest point in the image, and find the optimal orientation of any underlying lines  by radially integrating intensity.
        % 2) Once we have a start point and a direction, we continue searching for the underlying line by looking for intensity peaks until a threshold is met.
        % 3) We take this estimated microtubule, and the pixel region it spans and set it to the median value of the original image. We say that we have removed the intensity due to this microtubule.
        % 4) We repeat steps 1-3 until the microtubule estimates are smaller than a certain length
        % 5) We will fit a polynomial through each initial estimate (quadratic or
        %    cubic). This will be our final estimate for the microtubules.

        vis = estParams.visibility;
        fov = estParams.fieldOfView;
        step = estParams.stepSize;
        minMTlength = estParams.minMTlength;
        maxIter = estParams.maxMTsearchIterations;

        [~, ~, imSteer] = steerableDetector( imageGaussian, 4, 2);
        % remove very weak lines
        imSteerDet = imSteer.* imMask;
        steerTh = multithresh( imSteerDet(imMask), 2); imSteerDet( imSteerDet < steerTh(1) ) = 0;    
        imSteerDet = bwareafilt( logical(imSteerDet), [10, Inf]);

        intMedInit = median( imageGaussian( imMask) );
        TT = multithresh( imageGaussian(imMask), 1);
        imEstSteer = mat2gray( imageGaussian .* imSteerDet);
        imEst = imageGaussian.*imMask;
        dispImg( imEst, imSteer, imSteerDet, [ 1 3])
        intMax = 1.2*TT;
        % imEst = imageGaussian;
        iter = 0;
        coords = {};
        orients = [];
        maxIter=2;
        while intMax > TT && iter < maxIter 
             
            iter = iter +1
            % Find Bright Point in Image
            [ intMax, idxMax] = max( imEstSteer(:) );
            [ yMax, xMax] = ind2sub( size(imEst), idxMax);
            [coordTemp, success, orientStart] = estimateNextPoint( [ xMax; yMax], 0, step, vis, 540, imageGaussian);


            if success
                angOrient = [ orientStart, orientStart+pi];
                coordCurr = [];
                addedMT = 0;  
                for jAng = 1:2
                ang = angOrient( jAng);
                % Continue finding next points
                coordMT{jAng} = estimateMicrotubulePoints( [xMax; yMax], ang, imageGaussian, step, vis, fov);
                % Only save if mt length exceeds minimum acceptable estimate length
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
                imMT = imdilate( imMT, strel('disk', 3) );
                imEst( imMT) = intMedInit;
                minV = min( imEst(:) ); maxV = max( imEst(:) ); 
                imEst = mat2gray( imgaussfilt( imEst, 0.5) ) * ( maxV-minV) + minV;
                imEst = imEst .* imMask;
                imEstSteerOld = imEstSteer;
                imEstSteer( imMT) = 0;
                end
                % de-intensity the mtoc point
                spread=4;
                for jX = 1 : size(imEst, 2); for jY = 1 : size(imEst, 1)
                    if norm( [jX-xMax, jY-yMax] ) <= spread
                        imEst(jY, jX) = intMedInit;
                        imEstSteer(jY, jX)=0;
                    end
                end; end
                im2 = imEst;
            
            else
                imEstSteerOld = imEstSteer; 
                im1=imEst;
                spread=6;
                for jX = 1 : size(imEst, 2); for jY = 1 : size(imEst, 1)
                    if norm( [jX-xMax, jY-yMax] ) <= spread
                        imEst(jY, jX) = intMedInit;
                        imEstSteer(jY, jX)=0;
                    end
                end; end
                im2=imEst;
            end

            figure('NumberTitle', 'off', 'Name', ['Iter_' num2str(iter)] )
            subplot(231)
            imagesc(imageGaussian); colormap gray
            title('Original Image')
            subplot(233)
            imagesc( imEstSteerOld); colormap gray; hold on; plot( xMax, yMax, 'b*', 'MarkerSize', 14); hold off;
            title( 'Possible MTOCs' )
            subplot(232)
            imagesc(im1); colormap gray
            title( 'Image that is searched' )
            subplot(234)
            imagesc( im1); colormap gray; hold on; plot( coordMT{1}(1, :), coordMT{1}(2,:), 'r-', 'LineWidth', 2); plot( coordMT{2}(1, :), coordMT{2}(2,:), 'r-', 'LineWidth', 2); plot( xMax, yMax, 'b*', 'MarkerSize', 14); hold off
            title( 'Microtubule found' )
            subplot(235)
            imagesc(im2); colormap gray;
            title( 'Next Iter: Possible MT image')
            subplot(236)
            imagesc(imEstSteer); colormap gray
            title( 'Next Iter: Possible MTOC image')


        end

end

