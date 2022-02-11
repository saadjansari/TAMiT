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
                
                case 'Sid4'
                    feature = MitoticCell.findFeaturesDeNovo_Sid4( image, obj.params.estimate.sid4, obj.featureProps.spbBank);
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
            % Params
            params_spindle_finder.linewidth = 3;
            params_spindle_finder.brightestPixelAsSPB = 0;
            params_spindle_finder.spindleDeterminationSensitivity = 0.4;
            params_spindle_finder.spindleMinIntensity = 0.4;
            params_spindle_finder.visuals = 0;
            params_spindle_finder.visuals_path = '';
            params_spindle_finder.verbose = 0;
            spindle = MitoticCell.findTheSpindle( imageIn, params_spindle_finder);

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
        function spindle  = findTheSpindle( imageIn, params)
            % Find's a bright 2D spindle line
            
            imageIn_old = imageIn;
            
            % Params
            % Thickness of spindle axis used to calculate spindle intensity
            linewidth = params.linewidth;
            
            % The minimum normalized intensity that is used to find the end-points of
            % the spindle. The end-point is where the intensity function drops to this
            % value.
            spindleMinIntensity = params.spindleMinIntensity;
            
            % A setting (either 0 or 1) that forces one of the ends of the spindle to
            % be at the brightest pixel.
            brightestPixelAsSPB = params.brightestPixelAsSPB;
            
            % A sensitivity parameter used in the extended maxima function call. 
            spindleDeterminationSensitivity = params.spindleDeterminationSensitivity;
            
            % Other parameters
            visuals = params.visuals; % to turn on visuals
            visuals_path = params.visuals_path; % location to save the visuals
            verbose = params.verbose; % descriptions ON
            
            % Create visuals path if required
            if visuals == 1
                mkdir(visuals_path);
            end

            % Useful variables  
            % zAnisotropy = ceil( sizeVoxelsZ / sizeVoxelsX );
            numVoxelsX = size( imageIn, 1); numVoxelsY = size( imageIn, 2); numVoxelsZ = size( imageIn, 3);
            
            % Mask
            imageIn = mat2gray(imageIn);
            imMask3D = imageIn > 0;
            
            % Set background to min value (then rescale image to stretch
            % it)
            minVal = min( imageIn(imageIn(:) ~= 0) );
            imageIn = mat2gray((imageIn-minVal).*imMask3D);

            % Find Strong Signal Regions{{{
            % Convolve the image with a 3D gaussian to bring out the signal from the SPB
            image3DConv = imgaussfilt3( imageIn, 1, 'FilterDomain', 'spatial') .* imMask3D;
            imPlane = mat2gray( max( image3DConv, [], 3) );

            % Keep the strongest signal pixels
            % XXX (otsu thresholding doe not give correct separation of
            % desired voxels)
            % shift to multithresh (multiple otsu thresholds), pick first
            % threshold above the median value.
            threshOtsu = multithresh( imPlane( imPlane > 0), 1 );
            if threshOtsu(1) > median( imPlane( imPlane > 0) )
                threshOtsu = threshOtsu(1);
%             elseif threshOtsu(2) > median( imPlane( imPlane > 0) )
%                 threshOtsu = threshOtsu(2);
            else
                threshOtsu = median( imPlane( imPlane > 0) );
            end
            spindleMinIntensity = max( [threshOtsu spindleMinIntensity]);
            if verbose == 1
                 fprintf('Spindle Threshold / Spindle Min Intensity Value = %.3f\n',spindleMinIntensity);
            end

%             threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.4]);
            imPlaneStrong = imPlane;
            imPlaneStrong( imPlaneStrong < threshOtsu) = 0;

            imMask = imextendedmax( imPlaneStrong, spindleDeterminationSensitivity);
%             imPlaneStrongBool = imPlaneStrong;
%             imPlaneStrongBool( imPlaneStrongBool > 0 ) = 1;

            % }}}

            % Pick Brightest Strong Region and find its shape {{{

            % we'll pick the maxima which has the biggest mean intensity
            SS = regionprops( imMask, imPlane, 'Centroid', 'MajorAxisLength', 'Orientation', 'MeanIntensity','Area');
            SScell = struct2cell( SS);
            
            % Keep maxima with largest integrated intensity
            sumInts = cell2mat( SScell(1,:) ).*cell2mat( SScell(5,:) );
            [ spindleInt, idx] = max( sumInts);
            spindleInt = spindleInt/SScell{1,idx};

            if verbose == 1
                 fprintf( 'Spindle Mean Intensity = %.2f\n', spindleInt) 
            end
            if spindleInt < 0.6
                warning( 'The mean spindle intensity is less than 60% of the maximum intensity in the image. There may be issues with detection.');
            end

            % Now we can obtain the shape parameters of the spindle
            centroid = cell2mat( SScell( 2, idx) );
            majAxisLen = cell2mat(  SScell( 3, idx) );
            angleOrient = deg2rad( cell2mat(  SScell( 4, idx) ) );

            if angleOrient < 0
                angle = 3*pi/2 - abs(angleOrient);
            elseif angleOrient >=0
                angle = pi/2 + angleOrient;
            end

            if verbose == 1
                 fprintf( 'AngleRP = %.2f , AngleNew = %.2f\n', angleOrient, angle) 
            end

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

            % Measure Intensity
            IntSpindle = zeros(1, length(coords) );
            idxMaxPix = zeros(1, length(coords) );
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

            if visuals == 1

                % Figure 1
                h = figure; 
                imMaskRGB = zeros( numVoxelsX, numVoxelsY, 3); imMaskRGB(:,:,3) = imPlane; imMaskRGB(:,:,2) = imPlane.*imMask; imMaskRGB(:,:,1) = imPlane.*imMask;
                subplot(131); imagesc(imPlane); title('Original Image'); axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                subplot(132); imagesc(imPlaneStrong); title('Threshold (Otsu)');axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                subplot(133); imagesc(imMaskRGB); title('Extended Maximum');axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                %set(findall(gcf,'-property','FontSize'),'FontSize',16);
                saveas(h, [visuals_path,filesep,'preprocessing_step.png'])
                close(h)

                % Figure 2
                h = figure;
                subplot(131); imagesc(imPlane); title('Original Image'); axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                subplot(132); imagesc(imPlane); axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                hold on, line( [ ptMin(2), ptMax(2)], [ ptMin(1), ptMax(1)] , 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 10, 'Color', 'r' ); hold off; 
                title('Line through maximum'); 
                subplot(133), hold on;
                plot( IntSpindle, 'r-', 'LineWidth', 3); 
                yline( spindleMinIntensity, '-','Threshold','Color', 'k', 'LineWidth', 1);
                xl1 = xline( ind1, ':','Pole 1','Color', 'k', 'LineWidth', 2); 
                xl1.LabelVerticalAlignment = 'top';
                xl1.LabelHorizontalAlignment = 'center';
                xl1.LabelOrientation='horizontal';
                xl2 = xline( ind2, ':','Pole 2','Color', 'k', 'LineWidth', 2); 
                xl2.LabelVerticalAlignment = 'top';
                xl2.LabelHorizontalAlignment = 'center';
                xl2.LabelOrientation='horizontal';
                xlabel('Pixel number'); ylabel('Max intensity'); 
                title( sprintf('LineScan intensity with width %.1f', linewidth) ); hold off;
                %set(findall(gcf,'-property','FontSize'),'FontSize',16)
                saveas(h, [visuals_path,filesep,'estimation_step.png'])
                close(h)
                
                % Figure 3
                h = figure;
                imagesc(imPlane); colormap gray; axis equal;
                xlim([0,size(imageIn,2)]); ylim([0,size(imageIn,1)]);
                hold on, line( AstersX, AstersY, 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 8, 'Color', 'r' ), hold off;
                title('Final Estimated Spindle')
                %set(findall(gcf,'-property','FontSize'),'FontSize',16)
                saveas(h, [visuals_path,filesep,'estimation_final.png'])
                close(h)
            end

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
        
        % Sid4 {{{

        % findFeaturesDeNovo_Sid4 {{{
        function spbBank = findFeaturesDeNovo_Sid4( image2Find, params, props)
            % Mitotic Cell: Find Spindle Pole Bodies
            
            imageIn = im2double( image2Find);
            
            % Find all possible spots
            [ spb, intBkg] = MitoticCell.findSPB_sid4( imageIn, params);
            
            % Only keep the best 2 spots
            if length(spb) == 0 
                fprintf('WARNING!!!!!!!! No SPB spots found for this image!')
                fprintf('Returning empty spb_bank! User must ensure proper functioning!')
                spb_bank = {};
            end
            
            % Sort spbs in order of brightness
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
        % findSPB_sid4 {{{
        function [spbs, intBkg] = findSPB_sid4( imageOrg, params) 
            % Find SPB that are sid4-labeled

            % Params
%             visuals = 0;
%             visuals_path = '.';
%             verbose = 0;
%             min_spot_area = 25;
%             min_conn_region_area = 7;
            visuals = params.visuals;
            visuals_path = params.visuals_path;
            verbose = params.verbose;
            min_spot_area = params.min_spot_area;
            min_conn_region_area = params.min_conn_region_area;

            % Create visuals path if required
            if visuals == 1
                mkdir(visuals_path);
            end

            % Background
            intBkg = mean( imageOrg( imageOrg(:) > 0));
            
            % Get maximum Z projection and maxZ indices
            [im_max, idxZ] = max(imageOrg,[],3);
            
            % 3D mask
            imMaskB3 = logical( imageOrg);
         
            % gaussian filtered image (masked)
            imageG = mat2gray( imgaussfilt( imageOrg, 1) ); imageG = imageG .* imMaskB3;
            
            % Iterative thresholding
            imageG_max = max( imageG,[],3);
            tt = median( imageG_max( imageG_max>0) );
            accepted_area = 1000000; % initialized as SOME LARGE NUMBER 
            
            % Make movie
            if visuals == 1
                vidfile = VideoWriter( [visuals_path,filesep,'thresholding_movie.mp4'],'MPEG-4');
                vidfile.FrameRate = 20;
                open(vidfile);
                h = figure;
            end

            while tt < max(imageG_max(:)) && accepted_area >min_spot_area
                
                accepted_area = sum( imageG_max(:) > tt);
                tt = tt + 0.002;

                if verbose
                    disp(['threshold = ',num2str(tt),', N-true = ',num2str(accepted_area)]);
                end
                
                if visuals 
                    subplot(121)
                    imagesc(imageG_max); colormap gray; axis equal; 
                    xlim([0,size(imageG_max,2)]); ylim([0,size(imageG_max,1)]); 
                    title('Gauss filtered image')
                    subplot(122)
                    imagesc(imageG_max > tt); colormap gray; axis equal; 
                    xlim([0,size(imageG_max,2)]); ylim([0,size(imageG_max,1)]);
                    title(sprintf('Threshold = %.3f\n Num unmasked pixels = %d',tt,accepted_area))
                    writeVideo(vidfile, getframe(gcf));
                end
            end
            if visuals 
                close(vidfile)
                close(h)
            end
            
            img_threshed = bwareafilt( imageG_max>tt,[min_conn_region_area min_spot_area]);
            % We have extracted connected regions in the image of size
            % 9-20. Now, the goal is to find their centroids, and determine
            % if they are possible locations of sid4.
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
            
            if visuals == 1
                h = figure;
                subplot(121)
                imagesc(img_threshed); colormap gray; axis equal; 
                xlim([0,size(imageG_max,2)]); ylim([0,size(imageG_max,1)]); 
                title('Result after bwareafilt')
                subplot(122)
                imagesc(mat2gray(im_max)); colormap gray; axis equal; 
                xlim([0,size(imageG_max,2)]); ylim([0,size(imageG_max,1)]);
                hold on; plot(centroids(:,1), centroids(:,2), 'color','red', 'marker','o','linestyle','None', 'markersize',8, 'linewidth',2)
                title('Final spots position')
                saveas(h, [visuals_path,filesep,'final_estimate.png']);
                close(h)
            end
            % We might need to ensure that the spot found has an intensity
            % that is significantly bigger than the background intensity.
        end 
        % }}}
        % }}}
    end

end
