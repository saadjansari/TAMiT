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
        function obj = EstimateFeatures( obj, estimationImage, cTime, cChannel, idxChannel)
        % findFeatures : estimates and finds the features 
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Estimation
            obj.featureList{ idxChannel, cTime} = obj.EstimateFeaturesNovel( currentFeature, estimationImage);

        end
        % }}}

        % EstimateFeaturesNovel {{{
        function feature = EstimateFeaturesNovel( obj, currentFeature, image)

            switch currentFeature
                case 'Microtubule'
                    feature = MitoticCellBud.findFeaturesDeNovo_MT( image, obj.params.estimate.spindleBud, obj.featureProps.spindle_bud);
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

            % #################################
            % ###########  Params  ############
            % #################################

            displayFlag = params.display;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            spindleExclusionRange = params.spindleExclusionRange;
%             astralMTRepresentation = params.astralMTRepresentation;
%             astralMinLength = params.astralMinLength;
            sigma = params.common.sigma(1:dim);
            bkg = median( imageIn( imageIn(:) > 0) );
            spb_astral_dist = params.spb_astral_dist;
            voxel_size = params.voxel_size;
            curvature_model = params.curvature_model;
           
            % #################################
            % ###########  B K G   ############
            % #################################

            % Get background and nuclear background
            % Nuclear mask
            mask = MitoticCellBud.BY_find_nuclear_mask( imageIn);
            bkg_nuc = median( imageIn( find(mask(:) ) ) );

            % new background first
            bkg = median( imageIn( find(mask(:)==0 ) ) );

            % #################################
            % #########  Spindle Line #########
            % #################################

            % Set up parameters for spindle finder
            params_spindle_finder.linewidth = params.linewidth;
            params_spindle_finder.brightestPixelAsSPB = params.brightestPixelAsSPB;
            params_spindle_finder.spindleDeterminationSensitivity = params.spindleDeterminationSensitivity;
            params_spindle_finder.spindleMinIntensity = params.spindleMinIntensity;
            params_spindle_finder.visuals = params.visuals;
            params_spindle_finder.visuals_path = params.visuals_path;
            params_spindle_finder.verbose = params.common.verbose;

            % Extract Spindle. 
            spindle = MitoticCellBud.findTheSpindle( imageIn, params_spindle_finder);
            
            % Find amplitude along spindle line, then subtract background.
            lineAmpList = Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition );
            spindleAmp = median(lineAmpList) - bkg_nuc;
            
            % Create a Line feature for the Spindle MTs
            spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, spindleAmp, sigma, dim, props.fit{dim}.line, props.graphics.line);
            spindleMT.label = 'spindle';
%             spindleMT.SetBounds();

            % #################################
            % ########  Spindle Poles  ########
            % #################################
            
            AsterObjects{1} = {};
            AsterObjects{2} = {};
            
            if params.spindlePoles
                
                % Find amplitude at spindle pole locations
                spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) ) - spindleAmp;
                spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) ) - spindleAmp;
                
                % Force spindle pole body amplitude to be above atleast
                % equal to the background value
                if any(spbAmp < 0)
                    spbAmp( spbAmp < 0) = bkg;
                    if params.common.verbose > 1
                        warning( 'findFeaturesMT_deNovo : forcing SPBAmp to be > 0')
                    end
                end
                
                % Create the Spot features for SpindlePoleBodies
                SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); 
                SPB{1}.label = 'spb'; SPB{1}.SetBounds();
                SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); 
                SPB{2}.label = 'spb'; SPB{2}.SetBounds();
                AsterObjects{1} = {SPB{1}};
                AsterObjects{2} = {SPB{2}};
            end

            % #################################
            % ##########  Astral MTs  #########
            % #################################
            
            if params.astralMT

                % Find spindle angle w.r.t to each pole
                spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ) , 2*pi );
                spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

                % Get start position of astrals (spb_astral_dist away from the
                % spb opposite to the direction of spindle.
                spos = voxel_size.*spindleMT.startPosition; epos = voxel_size.*spindleMT.endPosition;
                dir = (epos-spos)/norm(epos-spos);
                sp{1} = spindleMT.startPosition+ (spb_astral_dist.*-dir)./voxel_size;
                sp{2} = spindleMT.endPosition+ (spb_astral_dist.*dir)./voxel_size;
                
                % Extract Astral Microtubules
                AMT{1} = MitoticCellBud.findAstralMicrotubules( imageIn, sp{1}, spindleAngle(1), spindleExclusionRange);
                AMT{2} = MitoticCellBud.findAstralMicrotubules( imageIn, sp{2}, spindleAngle(2), spindleExclusionRange); 

                % Create Curved MT objects
                for jAster = 1 : 2
                    idxRm = [];
                    curvedMTs = cell(1, length( AMT{jAster} ));
                    for jb = 1 : length( curvedMTs )
                        
                        coords = AMT{jAster}{jb};
                        
                        % get intial theta, curvature coefficients from
                        % coords
                        [origin,thetaInit,curvature, L, ccd] = get_tan_curvature_coeffs( coords, curvature_model);
                        
                        % Ensure origin is within the Z stacks
                        if origin(3) >= size( imageIn,3)
                            origin(3) = size(imageIn,3)-0.2;
                        elseif origin(3) <= 1
                            origin(3) = 1.2;
                        end
                        
                        % Get amplitude
                        A1 = smooth( findAmplitudeAlongCurveCoords( max(imageIn,[],3), round(ccd(1:2,:)) ) - bkg_nuc);
                        A1( A1 < 0) = 0; amp = median(A1);

                        % Ensure initial theta is pointed to inside the image region
                        if dim == 3
                            if origin(3) == 1
                                thetaInit(1,2) = pi/2 - 0.03;
                            elseif origin(3) == size(imageIn, 3)
                                thetaInit(1,2) = pi/2 + 0.03;
                            end
                        end
                        
                        % Ensure curvatures are reasonable to start with
                        % If sign is similar, ensure they are less than
                        % some overall value
%                         if sign(nV(1)) == sign(nV(2))
%                             if abs( nV(1) ) > 0.03 && abs( nV(2) ) > 0.0003
%                                 idxRm = [idxRm; jb];
%                             end
%                         end

                        % Create
                        curvedMTs{jb} = CurvedMT( curvature_model, origin', thetaInit, curvature, L, amp, sigma, dim, props.fit{dim}.curve, props.graphics.curve);
                    end
                    curvedMTs( idxRm) = [];
                    
                    if ~isempty(curvedMTs)
                        % ensure mts are within the image region
                        curvedMTs = restrict_curvedMT_inside_image( imageIn, curvedMTs);

                        % reduce mt length based on image intensity
                        curvedMTs = threhold_mt_length_using_image( imageIn, curvedMTs);
                        
                        idxRm = [];
                        for jb = 1 : length( curvedMTs )
                            if curvedMTs{jb}.L < 7
                                idxRm = [idxRm; jb];
                            end
                        end
                        curvedMTs( idxRm) = [];
                    
                        AsterObjects{jAster} = { AsterObjects{jAster}{:}, curvedMTs{:} };
                    end
                    
                end

            end
            
            % ##############################################
            % ##########  Setup Feature Heirarchy  #########
            % ##############################################
            
            % Aster MT Objects
            % SPBs + Curved MTs stored in an Aster object

            % Store basic objects in object Hierarchy
            % Aster objects stored in a class object: SPBs + Curved MTs stored in an Aster
            Asters = cell(1,2);
            for jAster = 1 : 2
                Asters{jAster} = Aster( dim, AsterObjects{jAster}{:} );
                Asters{jAster}.parameters = params;
            end

            % Spindle Feature
            spindleObj = SpindleNew( dim, imageIn, {spindleMT, Asters{:} }, props);
            spindleObj.findEnvironmentalConditions();
            
            % Try setting mask
            spindleObj.mask = logical( imageIn); 
            
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
                im_mip = max( img,[],3);
                
                % Mask
                mask = BY_find_nuclear_mask(frame);
                threshold_factor = 2;
                threshold = threshold_factor*median( img( find(mask(:) ) ) );

                %f = figure;
                %ha = tight_subplot( 1, length(curvedMTs), 0.01, 0.01,0.01);

                % For each curved mt, gets its coords
                for jf = 1 : length(curvedMTs)

                    % GetCoords
                    cc = curvedMTs{jf}.GetCoords();

                    % Get intensity at these coordiantes
                    int = Cell.findAmplitudeAlongCurveCoords( im_mip, round( cc(1:2,:)) );
                    tt = 1 : length(int);

                    % Find a mask for this specific feature by usign its simulated
                    % image, then binarizing it, dilating it.
%                     mask = curvedMTs{jf}.simulateFeature( size(frame) );
%                     mask = imbinarize( mask, 0.001);
%                     mask = imdilate( mask, strel('sphere',5) );
%                     %mask3 = repmat( mask, [1,1, size(frame,3)]);
%                     med = median( img( find( mask(:) ) ) );
%                     threshold = threshold_factor*med;

                    idxEnd = length(int);
                    if any( int > threshold)
                        if any( int < threshold)
                            idxEnd = 1+length(int) - find( flip(int)> threshold, 1, 'first');
                        end
                    else
                        idxEnd = 1;
                    end
                    curvedMTs{jf}.L = (idxEnd/length(int) ) *curvedMTs{jf}.L;
                end
            end
        end
        % }}}
        
        % findTheSpindle {{{
        function [spindle, h] = findTheSpindle( imageIn, params)
            % Find's a bright 2D spindle line
            % Input Parameters:
            % -----------------
            %   1. imageIn : 3D array (XYZ intensity data)
            %   2. params  : struct
            %       A. linewidth : float (Thickness of spindle axis used to
            %           calculate spindle intensity)
            %       B. spindleDeterminationSensitivity : float (Extended
            %           Region Sensitivity) (Range 0-1)
            %       C. spindleMinIntensity : float (Minimum Intensity of
            %           spindle) (Range 0-1)
            %       D. brightestPixelAsSPB: boolean (A setting that forces one 
            %           of the ends of the spindle to be at the brightest pixel.)
            %       E. visuals : boolean (turn on visuals)
            %       F. visuals_path : str (location to save the visuals)
            %       G. verbose : 0,1,2 (print relevant info to log)

            imageIn_old = imageIn;
            h=0;
            
            % Extract Params
            % Thickness of spindle axis used to calculate spindle intensity
            linewidth = params.linewidth;
            
            % A setting (either 0 or 1) that forces one of the ends of the spindle to
            % be at the brightest pixel.
            brightestPixelAsSPB = params.brightestPixelAsSPB;
            
            % Extended max sensitivity
%             spindleDeterminationSensitivity = params.spindleDeterminationSensitivity;
            
            % Minimum Intensity of spindle
            spindleMinIntensity = params.spindleMinIntensity;
            
            % Other parameters
            visuals = params.visuals; % to turn on visuals
            visuals_path = params.visuals_path; % location to save the visuals
            verbose = params.verbose; % descriptions ON
            
%             spindleDeterminationSensitivity = 0.4;
%             spindleMinIntensity = 0.65;
%             linewidth = 2;
%             brightestPixelAsSPB = 0;

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
            
            if verbose > 0
                 fprintf( ' - Spindle Mean Intensity = %.2f\n', spindleInt) 
            end
            if spindleInt < 0.6 && verbose>1
                warning( 'The mean spindle intensity is less than 60% of the maximum intensity in the image. There may be issues with detection.');
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
                error(' - Spindle found is not long enough. spindle min intensity is possibly too high. ');
            end
            ind1 = indSpindle(1); ind2 = indSpindle(end);
            
            if length(indSpindle) > 7
                ind1 = ind1+1;
                ind2 = ind2-1;
            end

            % Here X is the horizontal coordinate.
            [AstersY, AstersX] = ind2sub( size(imPlane), idxMaxPix( [ind1, ind2] ) );
            
            if AstersX(1) == AstersX(2) && AstersY(1) == AstersY(2)
                error('Only one point found for spindle')
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
                [~, ~, nms3(:,:,jZ), ~] = steerableDetector( imgaussfilt(imageIn(:,:,jZ),1), 4, 2);
            end
            
            % create a mask to apply to steerable image
            mask = imgaussfilt3(imageIn, 2);
            st = Methods.GetImageStats(mask,0);
            mask( mask < 1.0*st.Median) = 0; mask(mask ~=0) = 1;
            
            % Do a sweep radial sweep and look for peaks corresponding to
            % rods.
            im2 = max(imageIn,[],3);
            im2g = imgaussfilt(im2, 1);
            
            % Find threshold intensity for mts
            imVals = im2g( im2g ~= 0);
            minPkHeight = 1.0*median( imVals );
            mask2 = BY_find_nuclear_mask(imageIn);
            nbkg = median( imageIn( find(mask2(:) ) ) );
            
            % Seek Long lines
            [phiIntensity, phiValues] = Cell.radIntegrate2D( im2g, startPoint, 15, 30);
             warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi1 ] = findpeaks( phiIntensity, phiValues, ...
                'MinPeakHeight', 1.25*nbkg, 'MinPeakProminence', 0.25*nbkg);
            peakPhi1 = mod( peakPhi1, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )
            
            % remove angle belonging to spindle.
            spindleAngle = mod( spindleAngle, 2*pi);
            idxRm = find( abs(peakPhi1-spindleAngle) < spindleExclusionRange | ...
                abs(peakPhi1-spindleAngle+2*pi) < spindleExclusionRange | ...
                abs(peakPhi1-spindleAngle-2*pi) < spindleExclusionRange);
            peakPhi1( idxRm)=[];
            
            % Seek Short lines
            [phiIntensity, phiValues] = Cell.radIntegrate2D( im2g, startPoint, 5, 15);
            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi2 ] = findpeaks( phiIntensity, phiValues, ...
                'MinPeakHeight', 1.25*nbkg, 'MinPeakProminence', 0.25*nbkg);
            peakPhi2 = mod( peakPhi2, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )
            
            % remove angle belonging to spindle.
            spindleAngle = mod( spindleAngle, 2*pi);
            idxRm = find( abs(peakPhi2-spindleAngle) < spindleExclusionRange | ...
                abs(peakPhi2-spindleAngle+2*pi) < spindleExclusionRange | ...
                abs(peakPhi2-spindleAngle-2*pi) < spindleExclusionRange);
            peakPhi2( idxRm)=[];
            
            % combine unique long and short lines 
            peakPhi = concat_unique(peakPhi1,peakPhi2,0.2);
            
            % Default Values
            defaultVisibility = 10;
            defaultFieldOfView = 70;
            defaultStepSize = 4;
            minLength = 10;
            
            AstralMicrotubules = {};
            
            % Contrast editing
            %im22 = mean(imageIn,3);
            im_temp = imgaussfilt3(imageIn,2);
            med = median( im_temp(imageIn(:)~=0) );
            
%             nms3 = 0*imageIn;
%             im_tempp = 0*imageIn;
%             for jZ = 1 : size(imageIn,3)
%                 im_tempp(:,:,jZ) = imadjust(im_temp(:,:,jZ), [1*med, 5*med],[]);
%                 [~, ~, nms3(:,:,jZ), ~] = steerableDetector( im_tempp(:,:,jZ), 4, 4);
%             end
            
            for pp = 1 : length(peakPhi)
                
                imSearch = sum(nms3,3).*max(mask,[],3);
                % Find the max intensity start point in close neighborhood
                sx = startPoint(1); sy = startPoint(2); sp = [sx,sy];
                maxInt = imSearch( round(sy), round(sx));
                for jx = sx+[-1,0,1]
                    for jy = sy+[-1,0,1]
                        int = imSearch( round(jy), round(jx));
                        if int > maxInt
                            maxInt = int;
                            sp = [jx,jy];
                        end
                    end
                end

                % Set the two opposite angles for search
                cc = Methods.estimateCurveCoords( sp', peakPhi(pp), imSearch, defaultStepSize, defaultVisibility, defaultFieldOfView, 1);
                cc(1:2,1) = startPoint(1:2)';
                
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
        
        function mask = BY_find_nuclear_mask( img3)

            % Check if image is 3D or 2D
            if length( size(img3) ) == 3
                dim = 3;
            else
                dim=2;
            end

            if dim==3
                img = imgaussfilt3( img3, 1);
            else
                img = imgaussfilt( img3,1);
            end
            img = mat2gray(img);

            % Define function for contrast adjustsment
            func = 'histeq'; % imadjust,histeq or adapthisteq

            % HISTEQ
            imEdit = feval( func, img);

            %Threshold
            thresh = multithresh(imEdit(:),2);
            imBinary = imbinarize( imEdit, thresh(1) );

            % SE
            se3 = strel('disk', 3);
            im_filled = zeros( size(imBinary) );
            for jZ = 1 : size(img3,3)
                im_dilated = imdilate( imerode( imBinary(:,:,jZ), se3), se3);
                im_filled(:,:,jZ) = imfill( im_dilated, 'holes');
            end
            im_filled = mat2gray(im_filled);

            %montage({im2,imEdit,imBinary,im_dilated, im_filled, imColor},'Size',[1 6])

            mask = im_filled;
        end

        
    end

end
