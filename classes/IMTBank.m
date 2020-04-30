classdef IMTBank < OrganizerMaster
% A IMTBank is a higher level feature that sits inside an interphase cell. It is composed of Microtubule asters 
    properties
        polyOrder % polynomial order of curves
    end

    methods (Access = public)
        
        % IMTBank {{{
        function obj = IMTBank( dim, image, featureList, props2Fit)

            obj = obj@OrganizerMaster( dim, image, featureList, props2Fit, 'IMTBank');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj)

            % Get general environmental properties
            [vec, vecLabels] = getVecEnvironment( obj);

            % Loop over mt arrays and get their vectors. Remove the MTOCs from the fit vectors. 
            for jAster = 1 : obj.numFeatures
                [vec_mtarray, vecLabels_mtarray] = getVec( obj.featureList{ jAster}, obj.props2Fit.fit{obj.dim}.aster.curve);
                vec = [vec , vec_mtarray];
                vecLabels_mtarray = strcat( 'B', num2str( jAster), '_', vecLabels_mtarray);
                vecLabels = { vecLabels{:}, vecLabels_mtarray{:} };
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels);

            % Aster Vector
            % Create the vector for each subfeature MT_Array and send it to the objects for absorption
            for jb = 1 : obj.numFeatures

                strAster = [ 'B', num2str(jb), '_' ];

                % Take the vector and find indices associated with the mt_array. 
                idxA = find( ~cellfun( @isempty, strfind( vecLabels, strAster ) ) );
                vecA = vec( idxA );
                vecLabelsA = erase( vecLabels( idxA ), strAster );
                obj.featureList{ jb} = absorbVec( obj.featureList{ jb}, vecA, vecLabelsA );

            end
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            vecList = {};
            vecLabelsList = {};
            objList = {};

            % bundle vectors
            for jb = 1 : obj.numFeatures
                [ vecList{jb} , vecLabelsList{jb}] = getVec( obj.featureList{jb}, obj.props2Fit.fit{obj.dim}.aster.curve);
                objList{jb} = obj.featureList{jb};
            end
                
            % Prepend environmental parameters to each vec and vecLabel
            [vecE, vecLabelsE] = getVecEnvironment( obj, obj.props2Fit.fit{obj.dim}.Environment);
            
            for jFeat = 1 : length( vecList)
                vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
                vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
            end

        end
        % }}}
        
        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, Image2Find)
            % Searches for missing subfeatures and adds them to the featureList
            % We look for missing microtubule bundles            

            bkg = median( Image2Find( Image2Find(:) > 0) );
            minLength = 20;
            [Image2Find2D, idxZ] = max(Image2Find, [],3);
            if obj.dim == 2
                sigma=[1.2 1.2];
            elseif obj.dim == 3
                sigma=[1.2 1.2 1.0];
            end
            % Remove intensity for curves
            mask = obj.getMaskWithoutFeatures( Image2Find);
            
            % Run Steerable filter to detect curves
            imSteer3 = 0*Image2Find;
            for jj = 1 : size(Image2Find,3)
                [~, imSteer3(:,:,jj)] = Methods.FilterImageForCurveDetection( Image2Find2D);
            end
            imSteerMask2 = max( imSteer3.*mask, [],3);
            
            % Get image stats
            stats = Methods.GetImageStats( obj.image, 0);
            
            % Find a curve in the 2D steerable image
            [coord, ~, out] = Methods.FindCurve2D( imSteerMask2, 'ThreshInt', stats.ThreshHigh);
            
            % Return if useful curve not found
            if ~out.success || out.Length < minLength 
                successAdd = 0;
                return
            end
            
            % Interpolate to get fine coordinates
            L = sum( sqrt( diff( coord(1,:)).^2 + diff( coord(2,:)).^2 ) );
            [cX,cY] = Methods.InterpolateCoords( coord(1,:), coord(2,:), round( L/ length(coord(1,:)) ) );
            % Find the idx of max intensity for the origin
            idxAll = sub2ind( size(Image2Find2D), round(cY), round(cX) );
            [amp, idxMax] = max( Image2Find2D( idxAll) ); 
            zz = idxZ( round(cY(idxMax)), round(cX(idxMax) ));
            coords{1} = [ cX(idxMax:end); cY(idxMax:end); zz*ones(1,length(cX)-idxMax+1)];
            coords{2} = [ cX(idxMax:-1:1); cY(idxMax:-1:1); zz*ones(1,idxMax)];
            
            % Get length
            L(1) = sum( sqrt( diff( coords{1}(1,:)).^2 + diff( coords{1}(2,:)).^2 + diff( coords{1}(3,:)).^2 ) );
            L(2) = sum( sqrt( diff( coords{2}(1,:)).^2 + diff( coords{2}(2,:)).^2 + diff( coords{2}(3,:)).^2 ) );
     
            % Get Coeff
            cf1 = Bundle.estimatePolyCoefficients( coords{1}, [3 3 1], linspace(0,L(1), size(coords{1},2)));
            cf2 = Bundle.estimatePolyCoefficients( coords{2}, [3 3 1], linspace(0,L(2), size(coords{2},2)));

            % Get coordinates from coeffs
%             t1 = linspace(0,L(1),size(coords{1},2));
%             t2 = linspace(0,L(2),size(coords{2},2));
%             x1 = polyval( cf1{1}, t1); y1 = polyval( cf1{2}, t1);
%             x2 = polyval( cf2{1}, t2); y2 = polyval( cf2{2}, t2);
                
            origin = [cf1{1}(end), cf1{2}(end),cf1{3}(end)];
            if origin(3) >= size( Image2Find,3)
                origin(3) = size(Image2Find,3)-0.2;
            elseif origin(3) <= 1
                origin(3) = 1.2;
            end
            % Get initial tangent vector
            tanInit{1} = [cf1{1}(3), cf1{2}(3), 0];
            tanInit{2} = [cf2{1}(3), cf2{2}(3), 0];

            % Angle of tanInit
            thetaInit = zeros(2,2);
            thetaInit(1,1) = atan2( tanInit{1}(2), tanInit{1}(1) );
            thetaInit(2,1) = atan2( tanInit{2}(2), tanInit{2}(1) );
            thetaInit(1,2) = pi/2; thetaInit(2,2) = pi/2;

            % Normal Magnitude Coefficients
            nV = zeros(2,2);
            nV(1,1) = 2*(cf1{1}(3)*cf1{2}(2) - cf1{1}(2)*cf1{2}(3));
            nV(1,2) = 6*(cf1{1}(3)*cf1{2}(1) - cf1{1}(1)*cf1{2}(3));
            nV(2,1) = 2*(cf2{1}(3)*cf2{2}(2) - cf2{1}(2)*cf2{2}(3));
            nV(2,2) = 6*(cf2{1}(3)*cf2{2}(1) - cf2{1}(1)*cf2{2}(3));

            % Get amplitude along each bundle
            A1 = smooth( Cell.findAmplitudeAlongCurveCoords( Image2Find2D, round(coords{1}(1:2,:)) ) - bkg);
            A2 = smooth( Cell.findAmplitudeAlongCurveCoords( Image2Find2D, round(coords{2}(1:2,:)) ) - bkg);
            A1( A1 < bkg) = bkg; A2( A2 < bkg) = bkg;
            % Threshold amplitude
            thr1 = multithresh( A1(:), 1); thr2 = multithresh( A2(:), 1);
            % Find parameters vals when amplitude is above threshold
            idx1 = find(A1 >thr1, 1, 'last'); idx2 = find(A2 >thr2, 1, 'last');
            LO(1) = sum( sqrt( diff( coords{1}(1,1:idx1)).^2 + diff( coords{1}(2,1:idx1)).^2 + diff( coords{1}(3,1:idx1)).^2 ) );
            LO(2) = sum( sqrt( diff( coords{2}(1,1:idx2)).^2 + diff( coords{2}(2,1:idx2)).^2 + diff( coords{2}(3,1:idx2)).^2 ) );
            LO_mu = mean( LO); LO_mu = max( [2.0, LO_mu]);
            L(1) = max( [L(1), 8]);L(2) = max( [L(2), 8]);

            % Find amplitude and amplitude enhancement factor
            amp = median( [ A1(A1 <= thr1); A2(A2 <= thr2)] );
            ef = median( [ A1(A1 > thr1); A2(A2 > thr2)] )/amp;
            if ef < 1.5
                disp('enhancement factor is less than 1. set to 1.6')
                ef = 1.6;
            end
            if ef > 4
                disp('enhancement factor is greater than 4. set to 3.9')
                ef = 3.9;
            end

            % Ensure coefficients are within the image region
            if obj.dim == 3
                if origin(3) == 1
                    thetaInit(1,2) = pi/2 - 0.03;
                    thetaInit(2,2) = pi/2 - 0.03;
                elseif origin(3) == size(obj.image, 3)
                    thetaInit(1,2) = pi/2 + 0.03;
                    thetaInit(2,2) = pi/2 + 0.03;
                end
            end
            thetaInit = [ thetaInit(1,:) , thetaInit(2,:) ];
            nV = [ nV(1,:) , nV(2,:) ];
                
            % Create 
            bundle = BundleNew( origin, thetaInit,nV, LO_mu, L, ef, amp, sigma, obj.dim, obj.featureList{1}.props2Fit, obj.featureList{1}.display);
                
            if bundle.GetLength() <= minLength
                successAdd = 0;
                return
            end

            % Add the bundle 
            obj.addFeatureToList( bundle);
            % Add feature ID and update the featureMap
            bundle.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            successAdd = 1;
                            
        end
        % }}}
        
        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % For an interphase array, we can only remove  interphase microtubules
            % For each aster we will find a possible microtubule to remove, then we will pick the best of the possibilities and remove it from our feature list.
            % Look at the residual under each astral microtubule
            % Remove the one with the biggest unit residual per length.
            % If there are no microtubules to remove
            successRemove = 1;        

            % If no features to remove, exit the function
            if obj.numFeatures == 0
                successRemove = 0; return
            end

            % Find prominence of features
            p = zeros(1, obj.numFeatures);
            for idx = 1 : obj.numFeatures
                
                % find Prominence
                p(idx) = obj.featureList{ idx }.amplitude * ...
                    obj.featureList{ idx }.GetLength() * ...
                    mean( obj.featureList{ idx }.sigma);

            end
            
            % Find least prominent worst feature
            [~, idxRm ] = min(p);

            % Remove the worst microtubule
            obj.removeFeatureFromList( idxRm);
                
        end
        % }}}
        
        % Update3DFrom2DSpecific {{{
        function obj = Update3DFrom2DSpecific(obj, obj2D)
            % Update background and numAsters

            obj.findEnvironmentalConditions();
            
            % For each aster, set the z value of the seed and the
            % z coefficient of the curves
            for jf = 1 : obj.numFeatures
                seed = obj.featureList{jF}.featureList{1};
                [~,seed.position(3)] = max( obj.image( seed.position(2), seed.position(1), :) );
                seed.sigma(3) = 1.0;
                for jc = 2 : obj.featureList{jf}.numFeatures
                    obj.featureList{jf}.featureList{jc}.startPosition(3) = seed.position(3);
                    obj.featureList{jf}.featureList{jc}.sigma(3) = 1.0;
                    obj.featureList{jf}.featureList{jc}.cZ = 0.0;
                end
            end
                
        end
        % }}}
        
        % GetStructBasicFeatures {{{
        function feats = GetStructBasicFeatures(obj)

            for jf = 1 : obj.numFeatures
                next = 1;
                feats(jf).type = obj.featureList{jf}.type;
                switch feats(jf).type
                    case 'Spot'
                        feats(next).position = obj.featureList{jf}.position;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        next = next+1;

                    case 'Line'
                        feats(next).startPosition = obj.featureList{jf}.startPosition;
                        feats(next).endPosition = obj.featureList{jf}.endPosition;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        feats(next).length = obj.featureList{jf}.GetLength();
                        feats(next).orientation = obj.featureList{jf}.GetOrientation();
                        next = next+1;

                    case 'Curve'
                        feats(next).startPosition = obj.featureList{jf}.startPosition;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        feats(next).length = obj.featureList{jf}.GetLength();
                        feats(next).orientation = obj.featureList{jf}.GetOrientation();
                        coords = obj.featureList{jf}.GetCoords(); 
                        feats(jf).endPosition = coords(:,end);
                        next = next+1;

                    case 'AsterMT'
                end
            end
        end
        % }}}

        % estimateFeatureNextFrame {{{
        function objNew = estimateFeatureNextFrame( obj, imgNext)
            % estimate the features for the next frame using its image
            
            % We will depolymerize the curves from both ends until the overlap zone, we call this the seed
            % We will improve the localization of the seed by looking for max int pixels near it
            % We polymerize the curve from the seed end in the given orientation using a curve finder
            % This is our new guess
            %
            im2 = max( imgNext,[],3);
            objNew = obj.copyDeep();
            objNew.image = imgNext;

            % Graphics set up
            nX = size( im2, 2); nY = size( im2, 1); intMin = min( im2(:) ); intMax = max( im2(:) );
            imageSettings = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';


            figure; subplot(231)
            imagesc( max( obj.image, [],3)); hold on; obj.displayFeature( gca); title( 'Old Feature; Old Frame'); eval( imageSettings)

            subplot(232);
            imagesc( im2); hold on; obj.displayFeature( gca); title( 'Old Feature; New Frame'); eval(imageSettings)

            % Depolymerize to get seeds
            for jc = 1 : obj.numFeatures

                % Get seed coords
                t = obj.featureList{jc}.T;
%                 seeds{jc}.coords = obj.featureList{ jc}.GetCoords( linspace( max([t(1)/3, 0.03]), min([t(2)+(1-t(2))/3, 0.97]) ) );
                seeds{jc}.coords = obj.featureList{ jc}.GetCoords( linspace( max([t(1), 0.03]), min([t(2), 0.97]) ) );
                
                % Get seed orientation
                phi = mean( atan2( diff( seeds{jc}.coords(2,:) ) , diff( seeds{jc}.coords(1,:) ) ) );
                seeds{jc}.phi = [phi, phi+pi];

            end
            subplot(233);
            imagesc( im2); hold on;
            for jc = 1 : length(seeds)
                line( seeds{jc}.coords(1,:), seeds{jc}.coords(2,:), obj.featureList{jc}.display{:})
            end
            title( 'Old Seeds'); eval( imageSettings)

            % Improve seeds
            r = -4:4;
            t = linspace(0,1); 
            for jc = 1 : length(seeds)
                xx = seeds{jc}.coords(1,:); yy = seeds{jc}.coords(2,:);

                % Find a perp region around the seed
                xnew = xx + r'*cos( seeds{jc}.phi(1)+pi/2); 
                ynew = yy + r'*sin( seeds{jc}.phi(1)+pi/2);
                idx = sub2ind( size(im2), round(ynew) ,round(xnew));

                % Find intensity of pixels in this region and find the max intensity pixel perp to the direction of the seed.
                int = im2( idx(:) ); int = reshape( int, size(idx) ); [~, idxm] = max( int, [],1);
                idxFinal = sub2ind( size(idx), idxm, 1:size(idx,2)); 
                xf = xnew( idxFinal); yf = ynew(idxFinal);
                
                % Smooth coords by fitting a linear line
                fitX = polyfit( t, xf, 1); fitY = polyfit( t, yf, 1); 
                seeds{jc}.coordsNew = [ polyval( fitX,t) ; polyval( fitY, t)]; 
                phi = mean( atan2( diff( seeds{jc}.coordsNew(2,:) ) , diff( seeds{jc}.coordsNew(1,:) ) ) );
                seeds{jc}.phiNew = [phi, phi+pi];
            end
            subplot(234);
            imagesc( im2); hold on;
            for jc = 1 : length(seeds)
                line( seeds{jc}.coordsNew(1,:), seeds{jc}.coordsNew(2,:), obj.featureList{jc}.display{:})
            end
            title( 'New Seeds'); eval( imageSettings)

            % Polymerize seeds
            imSteer3 = 0*imgNext;
            for jj = 1 : size(imgNext,3)
                [~, imSteer3(:,:,jj)] = Methods.FilterImageForCurveDetection( imgNext(:,:,jj));
            end
%             imSteer = max( imSteer3,[],3);
%             [imGauss, imSteer] = Methods.FilterImageForCurveDetection( im2);
            for jc = 1 : length(seeds)
                ccz = obj.featureList{jc}.GetCoords(); ccz = ccz(3,:);
                idx1 = floor( min(ccz));
                idx2 = ceil( max(ccz));
                imSteerBundle = max( imSteer3(:,:,idx1:idx2),[],3);
                c4mSeed{1} = Methods.estimateCurveCoords( seeds{jc}.coordsNew(:,end), seeds{jc}.phiNew(1), imSteerBundle, 4, 10, 60, 0);
                c4mSeed{2} = seeds{jc}.coordsNew; c4mSeed{2}(:,1) = [];c4mSeed{2}(:,end) = [];
                c4mSeed{3} = Methods.estimateCurveCoords( seeds{jc}.coordsNew(:,1), seeds{jc}.phiNew(2), imSteerBundle, 4, 10, 60, 0);
                
                h = figure;
                imagesc( imSteerBundle); axis equal; colormap gray; hold on;
                line( c4mSeed{1}(1,:), c4mSeed{1}(2,:), 'Color', 'red', 'LineWidth', 2)
                line( c4mSeed{2}(1,:), c4mSeed{2}(2,:), 'Color', 'green', 'LineWidth', 2)
                line( c4mSeed{3}(1,:), c4mSeed{3}(2,:), 'Color', 'blue', 'LineWidth', 2)
                close(h)
                
                % Treat as single line and get ordered end to end coordinates
                cEE = [ flip( c4mSeed{3}, 2), c4mSeed{2}, c4mSeed{1} ]; 
            
                % Fit a polynomial
                dParam = 0;
                for jk = 1 : size( cEE,2) -1
                    dParam = [ dParam, dParam(end) + sqrt( diff( cEE( 1, jk:jk+1)).^2 + diff(cEE(2,jk:jk+1)).^2)]; 
                end
                dParam = dParam/max(dParam);
                objNew.featureList{jc}.cX = polyfit( dParam, cEE(1,:), obj.featureList{jc}.order(1));
                objNew.featureList{jc}.cY = polyfit( dParam, cEE(2,:), obj.featureList{jc}.order(2));
                t = sort( [ dParam( size(c4mSeed{3},2)), dParam( size(c4mSeed{3},2) + size(c4mSeed{2},2)) ] );
                t(1) = max( [ t(1) 0.03]);
                t(1) = min( [ t(1) 0.97]);
                t(2) = max( [ t(2) 0.03]);
                t(2) = min( [ t(2) 0.97]);
                objNew.featureList{jc}.T = t;
                
                % add 3rd dim
                if obj.dim == 3
                    cZ(1:2) = [ 0 seeds{jc}.coords(3,round(end/2))];
                    if cZ(2) == 1
                        cZ(1) = 0.5;
                    elseif cZ(2) == size(imgNext, 3) 
                        cZ(1) = -0.5;
                    end
                    objNew.featureList{jc}.cZ = cZ;
                end
            end

            subplot(235);
            imagesc( im2); hold on;
            for jc = 1 : length(seeds)
                cc = objNew.featureList{jc}.GetCoords();
                line( cc(1,:), cc(2,:), 'Color', 'g', 'LineWidth', 2)
            end
            title( 'New Bundles'); eval( imageSettings);

            subplot(236);
            imagesc( im2); hold on;
            obj.displayFeature( gca);
            for jc = 1 : length(seeds)
                cc = objNew.featureList{jc}.GetCoords();
                line( cc(1,:), cc(2,:), 'Color', 'g', 'LineWidth', 2)
            end
            title( 'Old and New Bundles'); eval( imageSettings);

        end
        % }}}
        
        function imfMask = getMaskWithoutFeatures(obj, img)
           % Obtain a mask with the features removed
                           
            imfMask = ones( size(img));

            % Simulate features one at a time
            for jb = 1 : obj.numFeatures
                imf = obj.featureList{jb}.simulateFeature( size(img) );
                imfMask = imfMask .* (imf < 0.2*max(imf(:) ));
            end
                   
        end

    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'IMTBank')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            for jFeat = 1 : length( S.featureList)
                featureList{ jFeat} = Bundle.loadFromStruct( S.featureList{ jFeat} ); 
            end

            obj = IMTBank( S.dim, S.image, featureList, S.props2Fit);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap();

        end
        % }}}
        
        % findIntBank{{{
        function intBankObj = findIntBank( imageIn, params, props)
            % Finds an interphase MT bank in the given image

            % Get dimensionality and background
            dim = length( size(imageIn) );
            bkg = median( imageIn( imageIn(:) > 0) );

            % Define sigma for this bank
            if dim == 2
                sigma=[1.2 1.2];
            elseif dim == 3
                sigma=[1.2 1.2 1.0];
            end

            % Filter image 
            [imG, ~] = Methods.FilterImageForCurveDetection( imageIn);

            % Find the curves
            vars = Methods.struct2cellvars(params);
            coords = Methods.FindCurves( imG, 'Plot', 0, vars{:}); 
            nBundles = length( coords);
            
%             figure;
%             subplot(2,1,1)
%             imagesc( max(imageIn,[],3)); colormap gray; axis equal; xlim([0 150]); ylim([0 150]); xticks([]); yticks([]); 
            % Create interphase bundles
            for jb = 1 : nBundles

                % Get length
                L(1) = sum( sqrt( diff( coords{jb}{1}(1,:)).^2 + diff( coords{jb}{1}(2,:)).^2 + diff( coords{jb}{1}(3,:)).^2 ) );
                L(2) = sum( sqrt( diff( coords{jb}{2}(1,:)).^2 + diff( coords{jb}{2}(2,:)).^2 + diff( coords{jb}{2}(3,:)).^2 ) );
                [cX1,cY1,cZ1] = Methods.InterpolateCoords3( coords{jb}{1}(1,:), coords{jb}{1}(2,:), coords{jb}{1}(3,:), round( L(1)/ length( coords{jb}{1}(1,:) ) ) );
                [cX2,cY2,cZ2] = Methods.InterpolateCoords3( coords{jb}{2}(1,:), coords{jb}{2}(2,:), coords{jb}{2}(3,:), round( L(2)/ length( coords{jb}{2}(1,:) ) ) );

                % Get Coeff
                cf1 = Bundle.estimatePolyCoefficients( [cX1;cY1;cZ1], [3 3 1], linspace(0,L(1),length(cX1 )));
                cf2 = Bundle.estimatePolyCoefficients( [cX2;cY2;cZ2], [3 3 1], linspace(0,L(2), length(cY2)));
                
                % Get coordinates from coeffs
                t1 = linspace(0,L(1),length(cX1 ));
                t2 = linspace(0,L(2),length(cX2 ));
                x1 = polyval( cf1{1}, t1); y1 = polyval( cf1{2}, t1);
                x2 = polyval( cf2{1}, t2); y2 = polyval( cf2{2}, t2);
                
                origin = [cf1{1}(end), cf1{2}(end),cf1{3}(end)];
                if origin(3) >= size( imageIn,3)
                    origin(3) = size(imageIn,3)-0.2;
                elseif origin(3) <= 1
                    origin(3) = 1.2;
                end
                
                % Get initial tangent vector
                tanInit{1} = [cf1{1}(3), cf1{2}(3), 0];
                tanInit{2} = [cf2{1}(3), cf2{2}(3), 0];

                % Angle of tanInit
                thetaInit = zeros(2,2);
                thetaInit(1,1) = atan2( tanInit{1}(2), tanInit{1}(1) );
                thetaInit(2,1) = atan2( tanInit{2}(2), tanInit{2}(1) );
                thetaInit(1,2) = pi/2; thetaInit(2,2) = pi/2;

                % Normal Magnitude Coefficients
                nV = zeros(2,2);
                nV(1,1) = 2*(cf1{1}(3)*cf1{2}(2) - cf1{1}(2)*cf1{2}(3));
                nV(1,2) = 6*(cf1{1}(3)*cf1{2}(1) - cf1{1}(1)*cf1{2}(3));
                nV(2,1) = 2*(cf2{1}(3)*cf2{2}(2) - cf2{1}(2)*cf2{2}(3));
                nV(2,2) = 6*(cf2{1}(3)*cf2{2}(1) - cf2{1}(1)*cf2{2}(3));
                
                % Get amplitude along each bundle
                A1 = smooth( Cell.findAmplitudeAlongCurveCoords( max(imageIn,[],3), round([cX1;cY1]) ) - bkg);
                A2 = smooth( Cell.findAmplitudeAlongCurveCoords( max(imageIn,[],3), round([cX2;cY2]) ) - bkg);
                A1( A1 < bkg) = bkg; A2( A2 < bkg) = bkg;
                % Threshold amplitude
                thr1 = multithresh( A1(:), 1); thr2 = multithresh( A2(:), 1);
                % Find parameters vals when amplitude is above threshold
                idx1 = find(A1 >thr1, 1, 'last'); idx2 = find(A2 >thr2, 1, 'last');
                LO(1) = sum( sqrt( diff( cX1(1:idx1)).^2 + diff( cY1(1:idx1)).^2 + diff( cZ1(1:idx1)).^2 ) );
                LO(2) = sum( sqrt( diff( cX2(1:idx2)).^2 + diff( cY2(1:idx2)).^2 + diff( cZ2(1:idx2)).^2 ) );
                LO_mu = mean( LO); LO_mu = max( [2.0, LO_mu]);
                L(1) = max( [L(1), 8]);L(2) = max( [L(2), 8]);


                % Find amplitude and amplitude enhancement factor
                amp = median( [ A1(A1 <= thr1); A2(A2 <= thr2)] );
                ef = median( [ A1(A1 > thr1); A2(A2 > thr2)] )/amp;
                if ef < 1.5
                    disp('enhancement factor is less than 1. set to 1.6')
                    ef = 1.6;
                end
                if ef > 4
                    disp('enhancement factor is greater than 4. set to 3.9')
                    ef = 3.9;
                end
                    
                % Ensure coefficients are within the image region
                if dim == 3
                    if origin(3) == 1
                        thetaInit(1,2) = pi/2 - 0.03;
                        thetaInit(2,2) = pi/2 - 0.03;
                    elseif origin(3) == size(imageIn, 3)
                        thetaInit(1,2) = pi/2 + 0.03;
                        thetaInit(2,2) = pi/2 + 0.03;
                    end
                end
                thetaInit = [ thetaInit(1,:) , thetaInit(2,:) ];
                nV = [ nV(1,:) , nV(2,:) ];
                
                % Create 
                bundles{jb} = BundleNew( origin, thetaInit,nV, LO_mu, L, ef, amp, sigma, dim, props.fit{dim}.curve, props.graphics.curve);
                
%                 subplot(2,nBundles,2*(nBundles-1)+jb-1)
%                 imagesc( max(imageIn,[],3)); colormap gray; axis equal; xlim([0 150]); ylim([0 150]); xticks([]); yticks([]); hold on;
%                 plot(x1,y1, 'c*', 'markerSize', 4, 'linewidth', 4); plot(x2,y2, 'c*', 'markerSize', 4, 'linewidth', 4);
%                 bundles{jb}.displayFeature(gca);
%                 plot(cX1(1), cY1(1), 'c*', 'markerSize', 12, 'linewidth', 3); hold off
                
            end
            
            % Create the interphase aster bank.
            intBankObj = IMTBank( dim, imageIn, bundles, props);
%             intBankObj.polyOrder = params.PolyOrder;
            intBankObj.findEnvironmentalConditions();
            intBankObj.forceInsideMask();

        end
        % }}}

    end
end
