classdef MitoticCellBud < Cell
% This is a specialized cell of superclass Cell that is used for a mitotic cell of BUdding yeast.
    properties
        phase % cell phase (e.g. prophase, metaphase, anaphase) if available
    end

    methods
        
        % MitoticCell {{{
        function obj = MitoticCellBud( imageData, featuresInChannels, channelsToFit, params, varargin)
        % MitoticCell : constructor function for MitoticCell object     

            obj = obj@Cell( imageData, featuresInChannels, channelsToFit, params, varargin{:}, 'Type', 'Mitosis');

        end
        % }}}

        % EstimateFeatures {{{
        function obj = EstimateFeatures( obj, estimationImage, cTime, cChannel, idxChannel, timeReverse, newEstimate)
        % findFeatures : estimates and finds the features 
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Get Start time 
            lifetime = obj.imageData.GetLifetime;
            if timeReverse
                startTime = lifetime(2);
            else
                startTime = lifetime(1);
            end

            % Novel Estimation for first frame
            if cTime == startTime || newEstimate
                obj.featureList{ idxChannel, cTime} = obj.EstimateFeaturesNovel( currentFeature, estimationImage);
            % Propagate old feature for later frames
            else 
                obj = obj.PropagateOldFeature( idxChannel, cChannel, cTime, timeReverse);
            end

            % Special Tasks 
            % If cut7 channel, try to get spindle information from MT channel
            if strcmp( obj.featuresInChannels{ idxChannel}, 'Cut7')
                % Get spindle information from microtubule channel if possible
                % (MICROTUBULE CHANNEL MUST BE ANALYZED FIRST FOR THIS TO WORK)
                mtChannel = find( strcmp( obj.featuresInChannels, 'Microtubule' ) );
                if ~isempty( mtChannel)
                    obj.featureList{ idxChannel, cTime}.harvestSpindleInfo( obj.featureList{mtChannel, cTime} );
                end
            end
        end
        % }}}

        % EstimateFeaturesNovel {{{
        function feature = EstimateFeaturesNovel( obj, currentFeature, image)
            disp('- DeNovo') 

            switch currentFeature
                case 'Microtubule'
                    feature = MitoticCellBud.findFeaturesDeNovo_MT( image, obj.params.estimate.spindle, obj.featureProps.spindle_bud);

                case 'Kinetochore'
                    feature = MitoticCellBud.findFeaturesDeNovo_KC( image, obj.params.estimate.kcbank);

                case 'Cut7'
                    feature = MitoticCellBud.findFeaturesDeNovo_Cut7( image, obj.params.estimate.cut7dist);
            end

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function spindleObj = findFeaturesDeNovo_MT( imageIn, params, props)

            % Mitotic Cell:
            %   Find the Spindle microtubule
            %   From each spindle pole, find astral microtubules

            displayFlag = params.display;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            spindleExclusionRange = deg2rad(45);

            % Get background and nuclear background
            % Nuclear mask
            mask = BY_find_nuclear_mask( imageIn);
            bkg_nuc = median( imageIn( find(mask(:) ) ) );

            % new background first
            bkg = median( imageIn( find(mask(:)==0 ) ) );
            
            if dim==2
                sigma = [2.5,2.5];
            else
                sigma = [2.5,2.5,1.2];
            end

            % Spindle
            % Find the Spindle. 
            spindle = MitoticCellBud.findTheSpindle( imageIn);

            % Create the Spindle MT
            if params.spindleMT
                lineAmpList = Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition );
                spindleAmp = median(lineAmpList) - bkg_nuc;
                spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, spindleAmp, sigma, dim, props.fit{dim}.line, props.graphics.line);
                spindleMT.label = 'spindle'; spindleMT.SetBounds();
            else
                fprintf('Not estimating spindle MT\n')
            end

            % Spindle Poles
            AsterObjects{1} = {};
            AsterObjects{2} = {};
            if params.spindlePoles
                % Create the SpindlePoleBody objects
                spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) ) - spindleAmp;
                spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) ) - spindleAmp;
                if any(spbAmp < 0)
                    spbAmp( spbAmp < 0) = bkg;
                    warning( 'findFeaturesMT_deNovo : forcing SPBAmp to be > 0')
                end
                SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); 
                SPB{1}.label = 'spb'; SPB{1}.SetBounds();
                SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); 
                SPB{2}.label = 'spb'; SPB{2}.SetBounds();
                AsterObjects{1} = {SPB{1}};
                AsterObjects{2} = {SPB{2}};
            end

            % Astral MTs
            if params.astralMT

                % Find the angle of this spindle w.r.t to each pole
                spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ) , 2*pi );
                spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

                % Get start position of astrals (3-5 pixels away from the
                % spb opposite to the direction of spindle.
                rr = 3;
                sp{1} = spindleMT.startPosition + rr*[cos(spindleAngle(1)+pi), sin(spindleAngle(1)+pi), 0];
                sp{2} = spindleMT.endPosition + rr*[cos(spindleAngle(2)+pi), sin(spindleAngle(2)+pi), 0];
                
                % Find Astral Microtubules
                AMT{1} = MitoticCellBud.findAstralMicrotubules( imageIn, sp{1}, spindleAngle(1), spindleExclusionRange);
                AMT{2} = MitoticCellBud.findAstralMicrotubules( imageIn, sp{2}, spindleAngle(2), spindleExclusionRange); 

                % Create Curved MT objects
                for jAster = 1 : 2
                    idxRm = [];
                    curvedMTs = cell(1, length( AMT{jAster} ));
                    for jb = 1 : length( curvedMTs )
                        
                        coords = AMT{jAster}{jb};
                        
                        % Get length
                        L = sum( sqrt( diff( coords(1,:)).^2 + diff( coords(2,:)).^2 + diff( coords(3,:)).^2 ) );
                        nInt1 = round( L/ length( coords(1,:) ) );
                        [cX1,cY1,cZ1] = Methods.InterpolateCoords3( coords(1,:), coords(2,:), coords(3,:), nInt1 );
                        % Get Coeff
                        cf1 = Bundle.estimatePolyCoefficients( [cX1;cY1;cZ1], [3 3 1], linspace(0,L,length(cX1 )));
                        % Get coordinates from coeffs
                        t1 = linspace(0,L,length(cX1 ));
                        x1 = polyval( cf1{1}, t1); y1 = polyval( cf1{2}, t1);

                        % Get origin
%                         origin = [cf1{1}(end), cf1{2}(end),cf1{3}(end)];
                        origin = coords(:,1);
                        if origin(3) >= size( imageIn,3)
                            origin(3) = size(imageIn,3)-0.2;
                        elseif origin(3) <= 1
                            origin(3) = 1.2;
                        end
                        % Get initial tangent vector and theta vector
                        tanInit{1} = [cf1{1}(3), cf1{2}(3), cf1{3}(1)];
                        thetaInit = [atan2( tanInit{1}(2), tanInit{1}(1) ), pi/2];
                        % Normal Magnitude Coefficients
                        nV = [ 2*(cf1{1}(3)*cf1{2}(2) - cf1{1}(2)*cf1{2}(3)), ...
                            6*(cf1{1}(3)*cf1{2}(1) - cf1{1}(1)*cf1{2}(3))];

                        % Get amplitude along each bundle
                        A1 = smooth( Cell.findAmplitudeAlongCurveCoords( max(imageIn,[],3), round([cX1;cY1]) ) - bkg_nuc);
                        A1( A1 < 0) = 0; amp = median(A1);

                        % Ensure coefficients are within the image region
                        if dim == 3
                            if origin(3) == 1
                                thetaInit(1,2) = pi/2 - 0.03;
                            elseif origin(3) == size(imageIn, 3)
                                thetaInit(1,2) = pi/2 + 0.03;
                            end
                        end
                        if sign(nV(1)) == sign(nV(2)) && (abs(nV(1)) > 0.01 && abs(nV(2)) > 0.005)
                            idxRm = [idxRm; jb];
                        end

                        % Create
                        curvedMTs{jb} = CurvedMT( origin', thetaInit, nV, L, amp, sigma, dim, props.fit{dim}.curve, props.graphics.curve);
                    end
                    curvedMTs( idxRm) = [];
                    
                    if ~isempty(curvedMTs)
                        % ensure mts are within the image region
                        curvedMTs = restrict_curvedMT_inside_image( imageIn, curvedMTs);

                        % reduce mt length based on image intensity
                        curvedMTs = threhold_mt_length_using_image( imageIn, curvedMTs);
                        
                        AsterObjects{jAster} = { AsterObjects{jAster}{:}, curvedMTs{:} };
                    end

                end

            end

            % Store basic objects in object Hierarchy
            % Aster objects stored in a class object: SPBs + Curved MTs stored in an Aster
            Asters = cell(1,2);
            for jAster = 1 : 2
                Asters{jAster} = Aster( dim, AsterObjects{jAster}{:} );
            end

            % Spindle Feature
            spindleObj = SpindleNew( dim, imageIn, {spindleMT, Asters{:} }, props);
            spindleObj.findEnvironmentalConditions();
            
            function curvedMTs = restrict_curvedMT_inside_image( frame, curvedMTs)
        
                nvx = size(frame,2);
                nvy = size(frame,1);
                nvz = size(frame,3);
                % restrict mts within the image region
                for jf = 1 : length(curvedMTs)

                    outside = 0;
                    % GetCoords
                    cc = curvedMTs{jf}.GetCoords();

                    % Identify coords that exceed the size of the image, and remove them
                    rmIdx_x = find( cc(1,:) < 1 | cc(1,:) > nvx);
                    rmIdx_y = find( cc(2,:) < 1 | cc(2,:) > nvy);
                    rmIdx_z = find( cc(3,:) < 1 | cc(3,:) > nvz);
                    rmIdx = union( union( rmIdx_x , rmIdx_y ), rmIdx_z);
                    if ~isempty(rmIdx)
                        outside = 1;
                        disp('This line estimate exceeds the image region. Restricting it inside')
                    end
                    cnt = 0;
                    while outside && cnt < 20
                        cnt = cnt+1;
                        curvedMTs{jf}.L = 0.95*curvedMTs{jf}.L;
                        cc = curvedMTs{jf}.GetCoords();
                        rmIdx_x = find( cc(1,:) < 1 | cc(1,:) > nvx);
                        rmIdx_y = find( cc(2,:) < 1 | cc(2,:) > nvy);
                        rmIdx_z = find( cc(3,:) < 1 | cc(3,:) > nvz);
                        outside = ~isempty( union( union( rmIdx_x , rmIdx_y ), rmIdx_z) );
                    end

                end

            end
            
            function curvedMTs = threhold_mt_length_using_image( frame, curvedMTs)
                % restrict mts within the image region

                frame = mat2gray(frame);
                %mask = create_astral_mask( frame);
                img = imgaussfilt(frame,1);
                %med = median( img( :) );

                threshold_factor = 1.1;

                %f = figure;
                %ha = tight_subplot( 1, length(curvedMTs), 0.01, 0.01,0.01);

                % For each curved mt, gets its coords
                for jf = 1 : length(curvedMTs)

                    % GetCoords
                    cc = curvedMTs{jf}.GetCoords();

                    % Get intensity at these coordiantes
                    int = Cell.findAmplitudeAlongCurveCoords( img, round( cc) );
                    tt = 1 : length(int);

                    % Find a mask for this specific feature by usign its simulated
                    % image, then binarizing it, dilating it.
                    mask = curvedMTs{jf}.simulateFeature( size(frame) );
                    mask = imbinarize( mask, 0.005);
                    mask = imdilate( mask, strel('disk',8,6) );
                    %mask3 = repmat( mask, [1,1, size(frame,3)]);
                    med = median( img( find( mask(:) ) ) );
                    threshold = threshold_factor*med;

                    idxEnd = length(int);
                    if any( int < threshold)
                        idxEnd = 1+length(int) - find( flip(int)> threshold, 1, 'first');
                    end
                    curvedMTs{jf}.L = (idxEnd/length(int) ) *curvedMTs{jf}.L;
                end
            end
        end
        % }}}
        
        % findTheSpindle {{{
        function [spindle, ax] = findTheSpindle( imageIn, ax)
            
            if nargin < 2 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end

            % Params
            spindleDeterminationSensitivity = 0.4;
            spindleMinIntensity = 0.65;
            linewidth = 2;
            brightestPixelAsSPB = 0;
            imMask3D = imageIn > 0;

            % Useful variables  
            % zAnisotropy = ceil( sizeVoxelsZ / sizeVoxelsX );
            numVoxelsY = size( imageIn, 1); numVoxelsX = size( imageIn, 2); numVoxelsZ = size( imageIn, 3);

            % Find Strong Signal Regions{{{
            % Convolve the image with a 3D gaussian to bring out the signal from the SPB
            image3DConv = imgaussfilt3( imageIn, 1, 'FilterDomain', 'spatial') .* imMask3D;
            imPlane = mat2gray( max( image3DConv, [], 3) );

            % Keep the strongest signal pixels
            threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.4]);
            imPlaneStrong = imPlane;
            imPlaneStrong( imPlaneStrong < threshOtsu) = 0;

            %imMask = imextendedmax( imPlaneStrong, spindleDeterminationSensitivity);
            imPlaneStrongBool = imPlaneStrong;
            imPlaneStrongBool( imPlaneStrongBool > 0 ) = 1;
            imMask = imPlaneStrongBool;

            % }}}

            % Pick Brightest Strong Region and find its shape {{{

            % we'll pick the maxima which has the biggest mean intensity
            SS = regionprops( imMask, imPlane, 'Centroid', 'MajorAxisLength', 'Orientation', 'MeanIntensity');
            SScell = struct2cell( SS);
            meanInts = cell2mat( SScell(4,:) );
            [ spindleInt, idx] = max( meanInts);

            if spindleInt < 0.5
                error( 'The mean spindle intensity is less than 50% of the maximum intensity in the image');
            else
                fprintf('\t Success: Spindle Mean Intensity = %.2f\n', spindleInt)
            end

            % Now we can obtain the shape parameters of the spindle
            centroid = cell2mat( SScell( 1, idx) );
            majAxisLen = cell2mat(  SScell( 2, idx) );
            angleOrient = deg2rad( cell2mat(  SScell( 3, idx) ) );

            if angleOrient < 0
                angle = 3*pi/2 - abs(angleOrient);
            elseif angleOrient >=0
                angle = pi/2 + angleOrient;
            end
%             disp( sprintf( 'AngleRP = %.2f , AngleNew = %.2f', angleOrient, angle) )

            % }}}

            % Find coordinates of a long line passing through the region {{{

            % Once you have the centroid and the angle of the spindle,  get a slice at
            % that angle through the centroid. The expectation is that the spindle is
            % directed along this slice and we should be able to find 2 points
            % corresponding to the SPBs where the intensity starts to dramatically fall
            % off.

            if majAxisLen < 3
                majAxisLen = 3;
            end

            ptMin = round( flip(centroid) + (5*majAxisLen) * [ cos(angle), sin(angle)  ] );
            ptMax = round( flip(centroid) - (5*majAxisLen) * [ cos(angle), sin(angle)  ] );

            coords1 = round( linspace( ptMin(1), ptMax(1), ceil(12*majAxisLen) ) );
            coords2 = round( linspace( ptMin(2), ptMax(2), ceil(12*majAxisLen) ) );
            coords = [coords1 ; coords2];

            % remove any column in coords whose either entry is less than or equal
            % to 0. or greater than or equal to numVoxelsX or numVoxelsY
            rmCol = union( find( coords(1, :) < 1 | coords(1, :) > numVoxelsX), ...
                find( coords(2, :) < 1 | coords(2, :) > numVoxelsY) );
            coords( :, rmCol) = [];

            % }}}

            % Find Intensity on line with a non-zero width {{{

            perpMatrix = [-linewidth : linewidth] .* [cos( angle+pi/2 ) ; sin( angle+pi/2)];

            % measure intensity
            IntSpindle = zeros(1, length(coords) );
            idxMaxPix = zeros(1, length(coords) );
%             dispImg( imPlane); hold on; title('How linescan is done with a line width'); set(gca, 'FontSize', 16)
%             plot( [coords(2,1) coords(2,end)], [coords(1,1) coords(1,end)], 'r-', 'LineWidth', 1)
            for jPix = 1 : length(IntSpindle)
                coordsPerp = round( coords(:, jPix) +  perpMatrix);
                coordsPerp( 1, find( coordsPerp(1,:) > numVoxelsY) ) = numVoxelsY; coordsPerp( 1, find( coordsPerp(1,:) < 1) ) = 1;
                coordsPerp( 2, find( coordsPerp(2,:) > numVoxelsX) ) = numVoxelsX; coordsPerp( 2, find( coordsPerp(2,:) < 1) ) = 1;
%                 plot( coordsPerp(2,:), coordsPerp(1,:), 'c-', 'LineWidth', 1)
                idxPerp = sub2ind( [numVoxelsY, numVoxelsX], coordsPerp(1, :), coordsPerp(2, :) );
                [IntSpindle( jPix), idxMaxInt] = max( imPlane(idxPerp) );
                idxMaxPix( jPix) = idxPerp( idxMaxInt);
            end
            IntSpindle = mat2gray( smooth( IntSpindle, 4) );

            % }}}

            % Find Spindle endpoints {{{
            % Now lets obtain the two SPB points. 
            indSpindle = find( IntSpindle > spindleMinIntensity);
            if length( indSpindle) == 1
                error('spindle found is not long enough. spindle min intensity is possibly too high. ');
            end
            ind1 = indSpindle(1); ind2 = indSpindle(end);
            
            if length(indSpindle) > 7
                ind1 = ind1+1;
                ind2 = ind2-1;
            end

            % Here X is the horizontal coordinate.
            [AstersY, AstersX] = ind2sub( size(imPlane), idxMaxPix( [ind1, ind2] ) );
            
            if AstersX(1) == AstersX(2) && AstersY(1) == AstersY(2)
                error('only one point found for spindle')
            end

            if brightestPixelAsSPB

                % Now find the maximum intensity location in this image plane and set
                % this as one of the two SPBs
                [ ~, indMax] = max( imPlane(:) );
                [SPB1y, SPB1x] = ind2sub( size(imPlane), indMax);

                % We'll keep this max intensity SPB, and for the second we'll choose
                % one from the two we found by picking the one furthest from this max
                % position.
                dist1 = norm( [ SPB1y, SPB1x] - [ AstersY(1), AstersX(1)] );
                dist2 = norm( [ SPB1y, SPB1x] - [ AstersY(2), AstersX(2)] );

                if dist1 > dist2
                    SPB2y = AstersY(1);
                    SPB2x = AstersX(1);
                else
                    SPB2y = AstersY(2);
                    SPB2x = AstersX(2);
                end

                AstersX = [SPB1x, SPB2x];
                AstersY = [SPB1y, SPB2y];

            end

            % find Z coordinates of the spindle
            [~, zIdx] = max( image3DConv, [], 3);
            AstersZ = [ zIdx( round(AstersY(1)), round(AstersX(1)) ), zIdx( round(AstersY(2)), round(AstersX(2)) ) ];
            spindlePosition = [AstersX' , AstersY', AstersZ'];

            spindle.MT.startPosition = spindlePosition(1, :);
            spindle.MT.endPosition = spindlePosition(2, :);
            spindle.MT.amplitude = mean( Cell.findAmplitudeAlongLine( imageIn, spindlePosition(1, :), spindlePosition(2, :) ) );
            % }}}

            % Display Plots {{{
            if nargout > 1
                line( AstersX, AstersY, 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 4, 'Color', 'g' )
            end

            plotFlag_SpindleFinder = 0;
            if plotFlag_SpindleFinder

                % Figure 1
                imMaskRGB = zeros( numVoxelsY, numVoxelsX, 3); imMaskRGB(:,:,1) = imPlane; imMaskRGB(:,:,2) = imPlane.*imMask; 
                dispImg( imPlane, imPlaneStrong, imMaskRGB, [1 3]);
                subplot(131); title('Original Image')
                subplot(132); title('Otsu Background Removal')
                subplot(133); title('Extended Maximum')
                set(findall(gcf,'-property','FontSize'),'FontSize',16);

                % Figure 2
                dispImg( imPlane, imPlane, [1 3])
                subplot(131); title('Original Image')
                subplot(132); hold on, line( [ ptMin(2), ptMax(2)], [ ptMin(1), ptMax(1)] , 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 10, 'Color', 'r' ); hold off; title('Line through possible spindle'); 
                subplot(133), hold on; plot( IntSpindle, 'r-', 'LineWidth', 2); line( [1 length(IntSpindle)], [ spindleMinIntensity, spindleMinIntensity], 'Color', 'k', 'LineWidth', 1); xlabel('Pixel number'); ylabel('Max intensity'); title( sprintf('LineScan intensity with width %.1f', linewidth) ); hold off;
                set(findall(gcf,'-property','FontSize'),'FontSize',16)
                
                % Figure 3
                dispImg( imPlane)
                hold on, line( AstersX, AstersY, 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 8, 'Color', 'r' ), hold off;
                title('Final Estimated Spindle')
                set(findall(gcf,'-property','FontSize'),'FontSize',16)
            end
            % }}}

        end
        % }}}
        
        % findAstralMicrotubules {{{
        function [AstralMicrotubules, ax] = findAstralMicrotubules( imageIn, startPoint, spindleAngle, spindleExclusionRange, ax)
            % Finds lines directed radially away from the start point
            % Removes the line that matches the spindle angle

            if nargin < 5 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end
            
            nms3 = 0*imageIn;
            for jZ = 1 : size(imageIn,3)
                [~, ~, nms3(:,:,jZ), ~] = steerableDetector(imageIn(:,:,jZ), 4, 3);
            end
            
            % create a mask to apply to steerable image
            mask = imgaussfilt3(imageIn, 2);
            st = Methods.GetImageStats(mask,0);
            mask( mask < 1.0*st.Median) = 0; mask(mask ~=0) = 1;
            
            % Do a sweep radial sweep and look for peaks corresponding to
            % rods.
            
            im2 = max(imageIn,[],3);
            im2g = imgaussfilt(im2, 1);
            [phiIntensity, phiValues] = Cell.radIntegrate2D( im2g, startPoint, 5, 15);
            
            % find peaks in Phi Intensity
            imVals = im2g( im2g ~= 0);
            minPkHeight = 1.0*median( imVals );

            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, ...
                'MinPeakHeight', minPkHeight, 'MinPeakProminence', 0.25*median( imVals ));
            peakPhi = mod( peakPhi, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )
            
            % remove angle belonging to spindle.
            spindleAngle = mod( spindleAngle, 2*pi);
            idxRm = find( abs(peakPhi-spindleAngle) < spindleExclusionRange | ...
                abs(peakPhi-spindleAngle+2*pi) < spindleExclusionRange | ...
                abs(peakPhi-spindleAngle-2*pi) < spindleExclusionRange);
            peakPhi( idxRm)=[];
            
            % Default Values
            defaultVisibility = 10;
            defaultFieldOfView = 70;
            defaultStepSize = 4;
            minLength = 10;
            
            AstralMicrotubules = {};
            
            for pp = 1 : length(peakPhi)
                               
                % Set the two opposite angles for search
                cc = Methods.estimateCurveCoords( startPoint(1:2)', peakPhi(pp), sum(nms3,3).*max(mask,[],3), defaultStepSize, defaultVisibility, defaultFieldOfView, 1);
                
                % figure out z-coordinates
                [a2D, i2d] = max( imageIn, [], 3);
                c3 = zeros( 1, size(cc,2));
                for jj = 1 : length(c3)
                    rng = 1; alist = []; idlist = [];
                    for j1 = -rng : rng
                        for j2 = -rng : rng
                            
                            if round(cc(2,jj))+j1 > 1 && round(cc(2,jj))+j1 < size(mask,1) && ...
                                    round(cc(1,jj))+j2 > 1 && round(cc(1,jj))+j2 < size(mask,2)
                                alist = [ alist, a2D( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                                idlist = [ idlist, i2d( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                            end
                            
                        end
                    end
                    [~, idd] = max( alist); c3( jj) = idlist(idd);
                end
                
                c3(1: ceil(length(c3)/2) ) = c3(1);
                
                try; c3(end) = mean( c3(end-2:end-1));end
                
                c3(ceil(length(c3)/2) : end) = mean( c3(ceil(2*length(c3)/3):end) );
                pz = polyfit( linspace(0,1, length(c3)), c3, 1);
                zz = polyval( pz, linspace(0,1, length(c3)));
                zz(zz < 1) = 1.2; zz(zz > size(imageIn,3) ) = size(imageIn,3)-0.2;
                cc = [cc; zz];
                
                % Length of end-end curve
                if norm( cc(:,end)-cc(:,1)) > minLength
                    AstralMicrotubules = { AstralMicrotubules{:}, cc};
                end
                
            end
            
            % Remove any duplicates (also remove if angle matches the
            % spindle angle)
            phiz = zeros(1, length(AstralMicrotubules));
            for j1 = 1 : length( AstralMicrotubules)
                phiz(j1) = atan2( AstralMicrotubules{j1}(2,2)-AstralMicrotubules{j1}(2,1) , AstralMicrotubules{j1}(1,2)-AstralMicrotubules{j1}(1,1) );
            end
            
            % spindle angle first
            idxRm = find( abs(phiz-spindleAngle) < spindleExclusionRange | abs(phiz-spindleAngle+2*pi) < spindleExclusionRange | abs(phiz-spindleAngle-2*pi) < spindleExclusionRange);
            AstralMicrotubules(idxRm) = []; phiz(idxRm) = [];
            
            destroy = [];
            for p1 = 1:length(phiz)
                for p2 = 1:length(phiz)
                    if p1 ~=p2
                        if abs( phiz(p1) - phiz(p2)) < 0.2 || abs( phiz(p1) - phiz(p2) -2*pi) < 0.2 || abs( phiz(p1) - phiz(p2) + 2*pi) < 0.2
                            % figure out which one to destroy
                            if norm( AstralMicrotubules{p1}(:,end)-AstralMicrotubules{p1}(:,1)) > norm( AstralMicrotubules{p2}(:,end)-AstralMicrotubules{p2}(:,1))
                                destroy = [destroy, p2];
                            else
                                destroy = [destroy, p1];
                            end
                        end
                    end
                end
            end
            destroy = unique(destroy);
            AstralMicrotubules( destroy) = [];

        end
        % }}}
        
        % }}}

        % Kinetochores {{{

        % findFeaturesDeNovo_KC {{{
        function kcBank = findFeaturesDeNovo_KC( image2Find, displayFlag)
            % Mitotic Cell: Find Kinetochores
            
            props.KC = {'position', 'amplitude', 'sigma'};
            display.KC = {'Color', [0 0.8 0] , 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 1};
            
            [ kc, maskNuclear, intNuclear] = MitoticCell.findKinetochores_mis12( image2Find);

            % Create the Kinetochores
            disp(sprintf( '                Number of kinetochores = %d', length(kc.x) ) )
            if length(kc.x) > 6
                error('too many kinetochores found. maybe you should fix the spot finder')
            end
            if length(kc.x) == 0 
                error('too little kinetochores found. maybe you should fix the spot finder')
            end
            
            dim = length( size( image2Find) ); 
            if dim==2, sigma=[1.2 1.2]; elseif dim==3, sigma=[1.2 1.2 1.0]; end

            for jspot = 1 : length(kc.x)
                kcBank{ jspot} = Spot( [ kc.x(jspot), kc.y(jspot), kc.z(jspot)], kc.amplitude(jspot)-intNuclear, sigma, dim, props.KC, display.KC);
                kcBank{ jspot}.findVoxelsInsideMask( logical(image2Find) );
            end

            % Create a Kinetochore Bank for handling and storage
            kcBank = KinetochoreBank( image2Find, kcBank{:} );
            kcBank.findEnvironmentalConditions();
            kcBank.syncFeaturesWithMap();
            kcBank.maskNuclear = maskNuclear;
            kcBank.backgroundNuclear = intNuclear - kcBank.background;

            if displayFlag
                f = figure;
                ax = axes;
                imagesc( max( image2Find, [], 3) ); colormap gray; axis equal; hold on
                kcBank.displayFeature( ax);
            end

        end
        % }}}
        % findKinetochores_mis12 {{{
        function [kc, maskNuclear, intNuclear] = findKinetochores_mis12( image2Find) 
            % Find Kinetochores that are mis12-gf labeled

            % Process:
            %   1. Find the 3D nucleus ( to a high accuracy)
            %   2. Find Kinetochores inside the nucleus by looking for peaks

            imageOrg = mat2gray(image2Find);
            imMaskB3 = logical( imageOrg);
            imMaskB2 = imMaskB3(:,:,1);
            [ image2D, idxMax] = max( image2Find, [], 3);

            imageG = mat2gray( imgaussfilt( imageOrg, 0.5) ); imageG = imageG .* imMaskB3;
            imVals = imageG( imageG(:) > 0); T = multithresh( imVals, 2);
            imMaskN3 = imageG; imMaskN3( imMaskN3 < T(1) ) = 0; imMaskN3 = logical( imMaskN3);

            % For each slice, dilate the image to get connected components, then only keep the biggest obejct, and then apply a convex hull
            for jZ = 1 : size( imMaskN3, 3)
                imSlice = imMaskN3(:,:,jZ);
                imSlice = imdilate( imSlice, strel('disk', 1) );
                [labeledImage, numberOfBlobs] = bwlabel( imSlice);
                blobMeasurements = regionprops(labeledImage, 'area');
                allAreas = [blobMeasurements.Area];
                [sortedAreas, sortIndexes] = sort(allAreas, 'descend');
                biggestBlob = ismember(labeledImage, sortIndexes(1));
                imSlice = biggestBlob > 0;
                imMaskN3(:,:,jZ) = bwconvhull( imSlice);

                %% MAYBE DRAW ELLIPSE AROUND BIG OBJECT

            end

            % Iterative thresholding
            image2D = image2D .*  max( imMaskN3, [], 3);
            imVal = image2D( image2D(:) > 0);
            threshVal = min( imVal);
            unmaskReg = sum( image2D(:) > threshVal);
            unmaskThresh = 25; % threshold until 25 pixels remain
            while unmaskReg > unmaskThresh
                unmaskReg = sum( image2D(:) > threshVal);
                threshVal = threshVal+0.01;
            end
            imgThreshed = image2D; imgThreshed( imgThreshed < threshVal) = 0;
            imgThreshed = imgThreshed .* bwconvhull( logical(imgThreshed) );
            
            % Find maxima in the unmasked regions
            maxima = imregionalmax( imgThreshed);
            [kc.y, kc.x] = find( maxima );
            % for each (x,y) pair find the z-pixel that has the max intensity, and find the amplitude
            kc.z = idxMax( sub2ind( size(image2Find), kc.y, kc.x) );
            kc.amplitude = image2Find( sub2ind( size(image2Find), round(kc.y), round(kc.x), round(kc.z) ) );
            maskNuclear = imMaskN3;
            intNuclear = sum( imMaskN3(:) .* image2Find(:) ) / sum( imMaskN3(:) );

        end
        % }}}
        
        % }}}

        % Cut7 {{{
        function Cut7 = findFeaturesDeNovo_Cut7( image2Find, displayFlag)
            % Find cut7

            Cut7 = Cut7Distribution( image2Find, {} );

        end
        % }}}

    end

end
