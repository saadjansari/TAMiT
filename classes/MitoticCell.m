classdef MitoticCell < Cell
% This is a specialized cell of superclass Cell that is used for a mitotic cell.
    properties
        phase % cell phase (e.g. prophase, metaphase, anaphase) if available
    end

    methods
        
        % MitoticCell {{{
        function obj = MitoticCell( image, lifetime, species, featuresInChannels, channelsToFit, settings)
        % MitoticCell : constructor function for MitoticCell object     

            obj = obj@Cell( image, lifetime, species, featuresInChannels, channelsToFit, 'Mitosis', settings);

        end
        % }}}

        % findFeatures {{{
        function obj = findFeatures( obj, cTime, cChannel, idxChannel)
        % findFeatures : estimates and finds the features 
            
            % FIRST FRAME
            if cTime == obj.lifetime(1)
                disp('            - DeNovo') 

                switch obj.featuresInChannels{ idxChannel}

                    case 'Microtubule'

                        obj.featureList{ idxChannel, cTime} = MitoticCell.findFeaturesDeNovo_MT( obj.image(:,:,:,cTime, cChannel), obj.settings.flags.debug, obj.settings.flags.fitSpindleOnly, obj.settings.flags.removeSpindleSPB);

                    case 'Kinetochore'

                        obj.featureList{ idxChannel, cTime} = MitoticCell.findFeaturesDeNovo_KC( obj.image(:,:,:,cTime, cChannel), obj.settings.flags.debug);

                    case 'Cut7'

                        obj.featureList{ idxChannel, cTime} = MitoticCell.findFeaturesDeNovo_Cut7( obj.image(:,:,:,cTime, cChannel), obj.settings.flags.debug);

                        % Also try getting spindle information from microtubule channel if it exists (MICROTUBULE CHANNEL MUST BE ANALYZED FIRST FOR THIS TO WORK)
                        mtChannel = find( strcmp( obj.featuresInChannels, 'Microtubule' ) );
                        if ~isempty( mtChannel)
                            obj.featureList{ idxChannel, cTime}.harvestSpindleInfo( obj.featureList{mtChannel, cTime} );
                        end

                end

            else % NOT FIRST FRAME

                disp('            - Using fits from previous timestep') 

                % Could be a bad frame, so we iterate back until we find a featureList ~= NaN that we can make a copy of
                usePrevFrame = 1;
                badFrame = 1;
                while badFrame
                    try
                        badFrame = isnan( obj.featureList{ idxChannel, cTime-usePrevFrame} )
                        usePrevFrame = usePrevFrame+1;
                    catch 
                        badFrame = 0;
                        obj.featureList{ idxChannel, cTime} = obj.featureList{ idxChannel, cTime-usePrevFrame}.copyDeep();
                        obj.featureList{ idxChannel, cTime}.image = obj.image(:,:,:,cTime, cChannel);
                    end
                end

            end

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function spindleObj = findFeaturesDeNovo_MT( imageIn, displayFlag, spindleOnlyFlag, removeSpindleSPB)

            % Mitotic Cell:
            %   Find the Spindle microtubule
            %   From each spindle pole, find astral microtubules

            displayFlag = displayFlag;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            spindleExclusionRange = deg2rad(45);

            props.SPB = {'position', 'amplitude', 'sigma'};
            props.SpindleMT = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            props.AsterMT = {'endPosition', 'amplitude', 'sigma'};
            props.Spindle = {'background'};
            display.SPB = {'Color', [0.7 0 0.7] , 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2};
            display.SpindleMT= {'Color', [1 0 1] , 'LineWidth', 5};
            display.AsterMT = {'Color', [1 0 1] , 'LineWidth', 3};
            if dim==3, sigma=[1.2 1.2 1.0]; elseif dim==2, sigma=[1.2 1.2]; end
            bkg = median( imageIn( imageIn(:) > 0) );

            % Spindle
            
            % Find the Spindle. 
            if displayFlag
                f = figure;
                ax = axes;
                imagesc( max(imageIn, [], 3) ); colormap gray; axis equal; hold on
                [spindle, ax]  = MitoticCell.findTheSpindle( imageIn, ax);
            else
                spindle = MitoticCell.findTheSpindle( imageIn);
            end
            % Create the Spindle MTs
            lineAmpList = Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition );
%             spindleAmp = lineAmpList( round( length(lineAmpList)/2) );
             spindleAmp = lineAmpList( round( length(lineAmpList)/2) ) - bkg;
%             spindleAmp = min( Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition ) )-bkg;
            spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, spindleAmp, sigma, dim, props.SpindleMT, display.SpindleMT );
            spindleMT.findVoxelsInsideMask( logical(imageIn) );

            % SPB
            
            % Create the SpindlePoleBody objects
%             spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) );
%             spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) );
             spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) )-bkg-spindleAmp;
             spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) )-bkg-spindleAmp;
            if any(spbAmp < 0);
                spbAmp( spbAmp < 0) = bkg;
                warning( 'findFeaturesMT_deNovo : forcing SPBAmp to be > 0')
            end
            SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), sigma, dim, props.SPB, display.SPB);
            SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), sigma, dim, props.SPB, display.SPB);

            % Astral MT
            if ~spindleOnlyFlag

                % Find the angle of this spindle w.r.t to each pole
                spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

                % Find Astral Microtubules
                if displayFlag 
                    [ spindle.Aster{1}.MT, ax] = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.startPosition, spindleAngle(1), spindleExclusionRange, ax);
                    [ spindle.Aster{2}.MT, ax] = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.endPosition, spindleAngle(2), spindleExclusionRange, ax); 
                else
                    spindle.Aster{1}.MT = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.startPosition, spindleAngle(1), spindleExclusionRange);
                    spindle.Aster{2}.MT = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.endPosition, spindleAngle(2), spindleExclusionRange); 
                end

                % Create Astral MTs
                AstralMT = cell(1,2);
                for jAster = 1 : 2
                    for jmt = 1 : length( spindle.Aster{jAster}.MT )

                        lineAmp = median( Cell.findAmplitudeAlongLine( imageIn, spindle.Aster{jAster}.MT{jmt}.startPosition, spindle.Aster{jAster}.MT{jmt}.endPosition ) )-bkg;

                        newMT = Line( spindle.Aster{jAster}.MT{jmt}.startPosition, spindle.Aster{jAster}.MT{jmt}.endPosition, lineAmp, sigma, dim, props.AsterMT, display.AsterMT);

                        newMT.findVoxelsInsideMask( logical(imageIn) );

                        if isempty( AstralMT{jAster} )
                            AstralMT{jAster} = {newMT};
                        else
                            AstralMT{jAster} = { AstralMT{jAster}{:}, newMT};
                        end
                    end
                end

            elseif spindleOnlyFlag % no astral microtubules allowed
                AstralMT{1} = {};
                AstralMT{2} = {};
            end

            % Store basic objects in object Hierarchy
            % SPBs + Astral MTs stored in an AsterMT
            % SpindleMT + 2 Asters stored in a Spindle
            for jAster = 1 : 2
                if ~isempty( AstralMT{jAster} )
                    AsterObjects{jAster} = AsterMT( dim, SPB{jAster}, AstralMT{jAster}{:} );
                else
                    AsterObjects{jAster} = AsterMT( dim, SPB{jAster} );
                end

            end
           
            if removeSpindleSPB
                spindleObj = Spindle( dim, imageIn, {spindleMT}, props.Spindle);
            else
                spindleObj = Spindle( dim, imageIn, {spindleMT, AsterObjects{:} }, props.Spindle);
            end
            spindleObj.findEnvironmentalConditions();
            spindleObj.syncFeaturesWithMap();
%             spindleObj.fillParams();

            %if displayFlag
                %spindleObj.displayFeature();
            %end
            
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
            spindleDeterminationSensitivity = 0.6;
            spindleMinIntensity = 0.75;
            linewidth = 3;
            brightestPixelAsSPB = 0;
            imMask3D = imageIn > 0;

            % Useful variables  
            % zAnisotropy = ceil( sizeVoxelsZ / sizeVoxelsX );
            numVoxelsX = size( imageIn, 1); numVoxelsY = size( imageIn, 2); numVoxelsZ = size( imageIn, 3);

            % Find Strong Signal Regions{{{
            % Convolve the image with a 3D gaussian to bring out the signal from the SPB
            image3DConv = imgaussfilt3( imageIn, 1, 'FilterDomain', 'spatial') .* imMask3D;
            imPlane = mat2gray( max( image3DConv, [], 3) );

            % Keep the strongest signal pixels
            threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.7]);
            imPlaneStrong = imPlane;
            imPlaneStrong( imPlaneStrong < threshOtsu) = 0;

            imMask = imextendedmax( imPlaneStrong, spindleDeterminationSensitivity);
            imPlaneStrongBool = imPlaneStrong;
            imPlaneStrongBool( imPlaneStrongBool > 0 ) = 1;

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
%                 disp('Success! A spindle was found')
                disp( sprintf( '               Success: Spindle Mean Intensity = %.2f', spindleInt) )
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

            coords1 = round( linspace( ptMin(1), ptMax(1), ceil(6*majAxisLen) ) );
            coords2 = round( linspace( ptMin(2), ptMax(2), ceil(6*majAxisLen) ) );
            coords = [coords1 ; coords2];

            % remove any column in coords whose either entry is less than or equal
            % to 0. or greater than or equal to numVoxelsX or numVoxelsY
            rmCol = [];
            for jCol = 1 : length(coords)
                
                if any( find( coords( :, jCol) < 1) )
                    rmCol = [ rmCol, jCol];
                elseif coords(1, jCol) > numVoxelsY
                    rmCol = [ rmCol, jCol];
                elseif coords(2, jCol) > numVoxelsX
                    rmCol = [ rmCol, jCol];
                end
                
            end
            % remove bad columns
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
                ind1 = indD-1; ind2 = indD+1;
            else
                ind1 = indSpindle(1); ind2 = indSpindle(end); end

            % Here X is the horizontal coordinate.
            [AstersY, AstersX] = ind2sub( size(imPlane), idxMaxPix( [ind1, ind2] ) );

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
                imMaskRGB = zeros( numVoxelsX, numVoxelsY, 3); imMaskRGB(:,:,1) = imPlane; imMaskRGB(:,:,2) = imPlane.*imMask; 
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
            % Finds straight lines directed radially away from the start point
            % Removes the line that matches the spindle angle

            if nargin < 5 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end
            
            Image2Find = imgaussfilt( imageIn, 1);

            if nargout > 1 
                [AstralMicrotubules, ax] = Cell.find_manyLines( Image2Find, startPoint, ax);
            else
                AstralMicrotubules = Cell.find_manyLines( Image2Find, startPoint);
            end

            if isempty( AstralMicrotubules)
                nomt = 1;
            else nomt = 0; end

            if ~nomt
                % remove the line associated with the spindle
                spindleAngleRange = spindleExclusionRange;

                for jmt = 1 : length( AstralMicrotubules)
%                     AstralMicrotubules{jmt}
                    AstralMicrotubulesAngles( jmt) = AstralMicrotubules{jmt}.phi;
                end
                
                % Difference of spindle MT angle and the suspected angle should be close to either 0, 2pi or -2pi
                angleDiff = abs( AstralMicrotubulesAngles-spindleAngle);
                angleDiff2= abs( angleDiff - 2*pi);

                idxSpindle = find( angleDiff < spindleAngleRange | angleDiff2 < spindleAngleRange);

                if ~isempty(idxSpindle) && nargout > 1
                    % Remove the MT line from the current axes
                    children = get(gca, 'Children');
                    childrenRev = children( length(AstralMicrotubules):-1:1);
                    delete( childrenRev( idxSpindle) );
                end

                AstralMicrotubules( idxSpindle) = [];

            end

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
