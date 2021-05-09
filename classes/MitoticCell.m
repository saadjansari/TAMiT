classdef MitoticCell < Cell
% This is a specialized cell of superclass Cell that is used for a mitotic cell.
    properties
        phase % cell phase (e.g. prophase, metaphase, anaphase) if available
    end

    methods
        
        % MitoticCell {{{
        function obj = MitoticCell( imageData, featuresInChannels, channelsToFit, params, varargin)
        % MitoticCell : constructor function for MitoticCell object     

            obj = obj@Cell( imageData, featuresInChannels, channelsToFit, params, varargin{:}, 'Type', 'Mitosis');

        end
        % }}}

        % EstimateFeatures {{{
        function obj = EstimateFeatures( obj, estimationImage, cTime, cChannel, idxChannel, timeReverse, newEstimate)
        % findFeatures : estimates and finds the features 
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Get Image for estimation
%             estimationImage = obj.imageData.GetImage;
%             estimationImage = estimationImage( :,:,:,cTime, cChannel);

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
                    feature = MitoticCell.findFeaturesDeNovo_MT( image, obj.params.estimate.spindle, obj.featureProps.spindle);

                case 'Kinetochore'
                    feature = MitoticCell.findFeaturesDeNovo_KC( image, obj.params.estimate.kcbank);

                case 'Cut7'
                    feature = MitoticCell.findFeaturesDeNovo_Cut7( image, obj.params.estimate.cut7dist);
                
                case 'Sid1'
                    feature = MitoticCell.findFeaturesDeNovo_Sid1( image, obj.params.estimate.sid1, obj.featureProps.spbBank);
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

            if dim==3, sigma=[1.2 1.2 1.0]; elseif dim==2, sigma=[1.2 1.2]; end
            bkg = median( imageIn( imageIn(:) > 0) );

            % Spindle
            % Find the Spindle. 
            spindle = MitoticCell.findTheSpindle( imageIn);

            % Create the Spindle MT
            if params.spindleMT
                lineAmpList = Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition );
%                 spindleAmp = lineAmpList( round( length(lineAmpList)/2) ) - bkg;
                spindleAmp = median(lineAmpList) - bkg;
                spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, spindleAmp, sigma, dim, props.fit{dim}.line, props.graphics.line);
                spindleMT.label = 'spindle';
            else
                fprintf('Not estimating spindle MT\n')
            end

            AsterObjects{1} = {};
            AsterObjects{2} = {};

            % Spindle Poles
            if params.spindlePoles
                % Create the SpindlePoleBody objects
                spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) ) - bkg - spindleAmp;
                spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) )-bkg-spindleAmp;
                if any(spbAmp < 0);
                    spbAmp( spbAmp < 0) = bkg;
                    warning( 'findFeaturesMT_deNovo : forcing SPBAmp to be > 0')
                end
                SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); SPB{1}.label = 'spb';
                SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); SPB{2}.label = 'spb';
                AsterObjects{1} = {SPB{1}};
                AsterObjects{2} = {SPB{2}};
            end

            % Astral MTs
            if params.astralMT

                % Find the angle of this spindle w.r.t to each pole
                spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ) , 2*pi );
                spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

                % Find Astral Microtubules
                AMT{1} = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.startPosition, spindleAngle(1), spindleExclusionRange);
                AMT{2} = MitoticCell.findAstralMicrotubules( imageIn, spindleMT.endPosition, spindleAngle(2), spindleExclusionRange); 

                % Create Astral MTs
                AstralMT = cell(1,2);
                for jAster = 1 : 2
                    idxRm=[]; 
                    mts = cell(1, length( AMT{jAster} ));
                    for jmt = 1 : length( mts )

                        lineAmp = median( Cell.findAmplitudeAlongLine( imageIn, AMT{jAster}{jmt}.startPosition, AMT{jAster}{jmt}.endPosition ) )-bkg;
                        mts{jmt} = Line( AMT{jAster}{jmt}.startPosition, AMT{jAster}{jmt}.endPosition, lineAmp, sigma, dim, props.fit{dim}.line, props.graphics.aster.line);
                        mts{jmt}.repr = 'spherical'; mts{jmt}.label = 'astralFY';
                        if mts{jmt}.GetLength() > 4
                            idxRm = [idxRm, jmt];
                        end
                    end
                    mts( idxRm) = [];
                    
                    if ~isempty( mts )
                        AsterObjects{jAster} = { AsterObjects{jAster}{:}, mts{:} };
                    end
                end
            end

            % Store basic objects in object Hierarchy
            
            % Aster MT Objects
            % SPBs + Astral MTs stored in an AsterMT
            Asters = cell(1,2);
            for jAster = 1 : 2
                Asters{jAster} = AsterMT( dim, AsterObjects{jAster}{:} );
            end

            % Spindle Feature
            % SpindleMT + 2 Asters stored in a Spindle
%             if isempty( AsterObjects) 
%                 spindleObj = Spindle( dim, imageIn, {spindleMT}, props);
%             else
                spindleObj = Spindle( dim, imageIn, {spindleMT, Asters{:} }, props);
%             end

            if displayFlag
%                 fName = 'estimate_spindle';
%                 f = figure( 'NumberTitle', 'off', 'Name', fName); 
%                 ax = axes;
%                 imagesc( ax, max(imageIn, [], 3) ); colormap gray; axis equal; hold on;
%                 spindleObj.displayFeature(ax)
                %sName = [fitInfo.saveDirectory, filesep, sName];
                %export_fig( sName, '-png', '-nocrop', '-a1') 
            end
            spindleObj.findEnvironmentalConditions();
            
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
            spindleMinIntensity = 0.7;
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
            threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.4]);
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
                
                image2D = max( im2double(imageIn), [], 3);

                % radially integrate in phi
                [ phiIntensity, phiValues] = Cell.radIntegrate2D( image2D, startPoint,3, 15);
    %             figure; plot( phiValues, phiIntensity, 'r-');

                % find peaks in Phi Intensity
                imVals = image2D( image2D ~= 0);
                mtBkg = median( imVals );
                minPkHeight = mtBkg + std( imVals );

                warning('off', 'signal:findpeaks:largeMinPeakHeight' )
                [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, 'MinPeakHeight', minPkHeight );
                peakPhi = mod( peakPhi, 2*pi);
                warning('on', 'signal:findpeaks:largeMinPeakHeight' )

                Lines = {};
                % Find length of lines
                for jLine = 1 : length( peakPhi)
                    cLine = Cell.find_singleLine3D( imageIn, startPoint, peakPhi( jLine) );
                    if ~isempty( cLine)
                        Lines = { Lines{:}, cLine}; 
                    end
                end
                AstralMicrotubules = Lines;
            end

            if isempty( AstralMicrotubules)
                nomt = 1;
            else nomt = 0; end

            if ~isempty( spindleAngle) && ~isempty( spindleExclusionRange) && ~nomt
                % remove the line associated with the spindle
                spindleAngleRange = spindleExclusionRange;

                for jmt = 1 : length( AstralMicrotubules)
%                     AstralMicrotubules{jmt}
                    AstralMicrotubulesAngles( jmt) = AstralMicrotubules{jmt}.phi;
                end
                
                % Difference of spindle MT angle and the suspected angle should be close to either 0, 2pi or -2pi
                angleDiff = abs( AstralMicrotubulesAngles-spindleAngle);
                angleDiff2= abs( angleDiff - 2*pi);
                angleDiff3= abs( angleDiff - 2*pi);
                
                idxSpindle = find( angleDiff < spindleAngleRange | angleDiff2 < spindleAngleRange | angleDiff3 < spindleAngleRange);

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

            imageIn = im2double(image2Find);
            Cut7 = Cut7Distribution( imageIn, {} );
            imVals = imageIn( imageIn(:) > 0);
            Cut7.background = median( imVals);

        end
        % }}}
        
        % Sid1 {{{

        % findFeaturesDeNovo_KC {{{
        function spbBank = findFeaturesDeNovo_Sid1( image2Find, params, props)
            % Mitotic Cell: Find Spindle Pole Bodies
            
            imageIn = im2double( image2Find);
            
            % Find all possible spots
            [ spb, intBkg] = MitoticCell.findSPB_sid1( imageIn);
            
            % Only keep the best 2 spots
            if length(spb) == 0 
                fprintf('WARNING!!!!!!!! No SPB spots found for this image!')
                spb_bank = cell( length(spb) );
                
            end
            
            if length(spb) > 2
                [~,idxKeep] = sort( [spb.amplitude], 'descend');
                idxKeep(3:end) = [];
                spb = spb( idxKeep);
            end
            fprintf( '\tNumber of SPBs found = %d\n', length(spb) )
            
            % Create the SPBs
            dim = length( size( image2Find) ); 
            if dim==2, sigma=[1.2 1.2]; elseif dim==3, sigma=[1.2 1.2 1.0]; end
            spb_bank = cell( 1, length(spb) );
            for jspot = 1 : length(spb)
                spb_bank{ jspot} = Spot( [ spb(jspot).x, spb(jspot).y, spb(jspot).z], spb(jspot).amplitude-intBkg, ...
                    sigma, dim, props.fit{dim}.spot, props.graphics.spot);
                spb_bank{ jspot}.findVoxelsInsideMask( logical(image2Find) );
                spb_bank{ jspot}.label = 'spb';
            end

            % Create a Spot Bank for handling and storage
            spbBank = SPBBank( dim, imageIn, spb_bank, props );
            spbBank.findEnvironmentalConditions();

%             if params.display
%                 f = figure;
%                 ax = axes;
%                 imagesc( max( image2Find, [], 3) ); colormap gray; axis equal; hold on
%                 spbBank.displayFeature( ax);
%             end

        end
        % }}}
        % findSPB_sid1 {{{
        function [spbs, intBkg] = findSPB_sid1( imageOrg) 
            % Find SPB that are sid1-labeled

            % Process:
            %   1. Find the 3D nucleus (to a high accuracy)
            %   2. Find Kinetochores inside the nucleus by looking for peaks

            intBkg = mean( imageOrg( imageOrg(:) > 0));
            
            % Get maximum Z projection and maxZ indices
            [im_max, idxZ] = max(imageOrg,[],3);
            
            % 3D mask
            imMaskB3 = logical( imageOrg);
         
            % gaussian filtered image (masked)
            imageG = mat2gray( imgaussfilt( imageOrg, 1) ); imageG = imageG .* imMaskB3;
            
            % Iterative thresholding
            imageG_max = max( imageG,[],3);
            tt = min( imageG_max(:));
            accepted_area = 1000000;
            min_spot_area = 25;
            
            % Make movie
%             vidfile = VideoWriter('~/Desktop/thesholding_movie.mp4','MPEG-4');
%             vidfile.FrameRate = 10;
%             open(vidfile);
%             figure;
            while tt < max(imageG_max(:)) && accepted_area >min_spot_area
                
%                 imagesc( [imageG_max , imageG_max > tt])
%                 title(['threshold = ',num2str(tt),', N-true = ',num2str(accepted_area)])
%                 
                accepted_area = sum( imageG_max(:) > tt);
%                 disp(['threshold = ',num2str(tt),', N-true = ',num2str(accepted_area)]);
                tt = tt + 0.001;
                
%                 ff = getframe(gcf);
%                 writeVideo(vidfile, ff);
            end
%             close(vidfile)
            
            img_threshed = bwareafilt( imageG_max>tt,[7 min_spot_area]);
            % We have extracted connected regions in the image of size
            % 9-20. Now, the goal is to find their centroids, and determine
            % if they are possible locations of sid1.
            
            s  = regionprops( img_threshed, im_max, 'Centroid', 'Area', 'MeanIntensity');
            centroids = cat(1, s.Centroid);
            spot_intensities = cat(1, s.MeanIntensity);
            spbs = [];
            
            for jrow = 1 : length( spot_intensities)
                spbs(jrow).x = centroids(jrow,1);
                spbs(jrow).y = centroids(jrow,2);
                spbs(jrow).z = idxZ( round(centroids(jrow,2)), round(centroids(jrow,1)) );
                spbs(jrow).amplitude = spot_intensities( jrow);
            end
            
%             figure;
%             imagesc( [im_max,im_max]); colormap gray; axis equal
%             hold on; plot(centroids(1), centroids(2), 'color','green', 'marker','*','linewidth',1, 'markersize',12)
            % We might need to ensure that the spot found has an intensity
            % that is significantly bigger than the background intensity.
        end 
        % }}}
        % }}}
    end

end
