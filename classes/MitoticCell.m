classdef MitoticCell < Cell
% This is a specialized cell of superclass Cell that is used for a mitotic cell.
    properties
        phase % cell phase (e.g. prophase, metaphase, anaphase) if available
    end

    methods
        
        % MitoticCell {{{
        function obj = MitoticCell( image, lifetime, species, featuresInChannels, settings)
        % MitoticCell : constructor function for MitoticCell object     

            obj = obj@Cell( image, lifetime, species, featuresInChannels, 'Mitosis', settings);

        end
        % }}}

        % findFeatures {{{
        function obj = findFeatures( obj, image, cTime, cChannel)
        % findFeatures : estimates and finds the features 
            

            if cTime == obj.lifetime(1)
                disp('            - DeNovo') 

                switch obj.featuresInChannels{ cChannel}

                    case 'Microtubule'

                        obj.featureList{ cChannel, cTime} = MitoticCell.findFeaturesDeNovo_MT( image );

                    case 'Kinetochore'

                        obj.featureList{ cChannel, cTime} = MitoticCell.findFeaturesDeNovo_KC( image );

                end

            else

                disp('            - Using fits from previous timestep') 
                obj.featureList{ cChannel, cTime} = obj.featureList{ cChannel, cTime-1}.copyDeep( );
                obj.featureList{ cChannel, cTime} = obj.featureList{ cChannel, cTime-1}.fillParams();

            end

        end
        % }}}

        % prepareFit {{{
        function [ fitProblem, fitInfo ] = prepareFit( obj, featureMain, imageOrg, parameters, fitScope, cFeature)
            
            % Input Checks
            if length( featureMain) ~= 1
                error('prep_fitLocal : there should only be a single main feature'), end
            
            if ~strcmp( fitScope, 'local') && ~strcmp( fitScope, 'global') && ~strcmp( fitScope, 'globum_add') && ~strcmp( fitScope, 'globum_remove')

                error('prepareFit : input argument fitScope must be either ''local'' or ''global'' '); end

            if strcmp( fitScope, 'local') && nargin < 6
                error('prepareFit : for a local fit, feature number argument must be provided'); end

            % The fit will be 3-dimensional
            imageOrg = im2double( imageOrg); 
            
            % set general optimization options 
%             opts = optimoptions( @lsqnonlin, ...
%                 'MaxFunEvals', 2000, ...
%                 'OptimalityTolerance', 1e-12, ...
%                 'MaxIter', 10, ...
%                 'TolFun', 1e-7, ...
%                 'FiniteDifferenceStepSize', 1e-1, ...
%                 'FiniteDifferenceType', 'central', ...
%                 'StepTolerance', 1e-3, ...
%                 'display', 'iter', ... 
%                 'OutputFcn', @plotFit );

            opts = optimoptions( @lsqnonlin, ...
                'MaxIter', 15, ...
                'TolFun', 1e-7, ...
                'FiniteDifferenceStepSize', 1e-2, ...
                'display', 'off', ...
                'OutputFcn', @plotFit );

            % Prompt the main feature for a list of subfeatures, their fit vectors, and their label vectors
            if strcmp(fitScope, 'local')
                [ fitVec, fitLabels, fitObj] = getVecLocal( featureMain);
                fitVec = fitVec{ cFeature};
                fitLabels = fitLabels{ cFeature};
                fitObj = fitObj{ cFeature};
            elseif strcmp(fitScope, 'global')
                [ fitVec, fitLabels ] = getVec( featureMain); 
                fitObj = featureMain; cFeature = 1; 
            elseif strcmp(fitScope, 'globum_add')
                [ fitVec, fitLabels ] = getVec( featureMain); 
                fitObj = featureMain; cFeature = 1; 
                fitInfo.Nnew = featureMain.getSubFeatureNumber();
                fitInfo.Nold = fitInfo.Nnew - 1;
            elseif strcmp(fitScope, 'globum_remove')
                [ fitVec, fitLabels ] = getVec( featureMain); 
                fitObj = featureMain; cFeature = 1; 
                fitInfo.Nnew = featureMain.getSubFeatureNumber();
                fitInfo.Nold = fitInfo.Nnew + 1;
            end


            
            % find upper and lower bounds of parameters
            [fitVecUb, fitVecLb] = MitoticCell.getUpperLowerBounds( fitVec, fitLabels, imageOrg);

            % Exploration Speed
            speedVec = MitoticCell.getExplorationSpeedVector( fitLabels);
            fitVec = fitVec./speedVec;
            fitVecUb = fitVecUb./speedVec;
            fitVecLb = fitVecLb./speedVec;
%             error('stop here')

            % Set up parameters for fit
            fitInfo.featureMain = featureMain;
            fitInfo.featureCurrent = fitObj;
            fitInfo.vecInit = fitVec;
            fitInfo.labels = fitLabels;
            fitInfo.mask = logical( imageOrg);
            fitInfo.featureIndex = cFeature;
%             fitInfo.parentObj = parameters.Cell;
%             fitInfo.Image2Fit = uint8(imageOrg);
            fitInfo.numVoxels.X = size( imageOrg,2);
            fitInfo.numVoxels.Y = size( imageOrg,1);
            fitInfo.numVoxels.Z = size( imageOrg,3);
            fitInfo.speedVec = speedVec;
            fitInfo.channel = parameters.channel;
            fitInfo.time = parameters.time;
            fitInfo.fitScope = fitScope;
            
            % Find voxel indices for gaussian feature computation
            % simulate feature
%             imSim = fitObj.simulateFeature( 0*imageOrg);
%             idx = find( imSim > 0.001*max( imSim(:) ) );
%             [y x z] = ind2sub( size(imageOrg), idx);
%             fitInfo.fastComputation.idx = idx;
%             fitInfo.fastComputation.x = x;
%             fitInfo.fastComputation.y = y;
%             fitInfo.fastComputation.z = z;

            % make the error function for lsqnonlin
            f = Cell.makeErrorFcn( imageOrg, fitInfo );

            % Crate Optimization Problem
            fitProblem = createOptimProblem( 'lsqnonlin', 'objective', f, 'x0', fitVec, 'ub', fitVecUb, 'lb', fitVecLb, 'options', opts);

            % Create handle for plotting function
            function stop = plotFit( x, optimV, state)
                stop = Cell.plot_midFit(x, optimV, state, fitInfo);
            end

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function spindleObj = findFeaturesDeNovo_MT( imageIn)

            % Mitotic Cell:
            %   Find the Spindle microtubule
            %   From each spindle pole, find astral microtubules

            plotFlag_estimate = 0;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            spindleExclusionRange = deg2rad(45);

            props.SPB = {'position', 'amplitude', 'sigma'};
            props.SpindleMT = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            props.AsterMT = {'endPosition', 'amplitude', 'sigma'};
            props.Spindle = {'background'};
            display.SPB = {'Color', 'Red', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2};
            display.SpindleMT= {'Color', 'Red', 'LineWidth', 5};
            display.AsterMT = {'Color', 'Red', 'LineWidth', 3};
            if dim==3, sigma=[1.2 1.2 1.0]; elseif dim==2, sigma=[1.2 1.2]; end
            bkg = median( imageIn( imageIn(:) > 0) );

            % Spindle
            
            % Find the Spindle. 
            if plotFlag_estimate
                f = figure;
                ax = axes;
                imagesc( max(imageIn, [], 3) ); colormap gray; axis equal; hold on
                [spindle, ax]  = MitoticCell.findTheSpindle( imageIn, ax);
            else
                spindle = MitoticCell.findTheSpindle( imageIn);
            end
            % Create the Spindle MTs
            spindleAmp = median( Cell.findAmplitudeAlongLine( imageIn, spindle.MT.startPosition, spindle.MT.endPosition ) )-bkg;
            spindleMT = Line( spindle.MT.startPosition, spindle.MT.endPosition, spindleAmp, sigma, dim, imageIn, props.SpindleMT, display.SpindleMT );

            % SPB
            
            % Create the SpindlePoleBody objects
            spbAmp(1) = imageIn( spindleMT.startPosition(2), spindleMT.startPosition(1), spindleMT.startPosition(3) )-bkg-spindleAmp;
            spbAmp(2) = imageIn( spindleMT.endPosition(2), spindleMT.endPosition(1), spindleMT.endPosition(3) )-bkg-spindleAmp;
            SPB{1} = Spot( spindleMT.startPosition, spbAmp(1), sigma, dim, imageIn, props.SPB, display.SPB);
            SPB{2} = Spot( spindleMT.endPosition, spbAmp(2), sigma, dim, imageIn, props.SPB, display.SPB);

            % Astral MT
            
            % Find the angle of this spindle w.r.t to each pole
            spindleAngle(1) = mod( atan2( spindleMT.endPosition(2)-spindleMT.startPosition(2) , spindleMT.endPosition(1)-spindleMT.startPosition(1) ), 2*pi );
            spindleAngle(2) = mod( atan2( spindleMT.startPosition(2)-spindleMT.endPosition(2) , spindleMT.startPosition(1)-spindleMT.endPosition(1) ) , 2*pi );

            % Find Astral Microtubules
            if plotFlag_estimate
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
                    newMT = Line( spindle.Aster{jAster}.MT{jmt}.startPosition, spindle.Aster{jAster}.MT{jmt}.endPosition, lineAmp, sigma, dim, imageIn, props.AsterMT, display.AsterMT);
                    if isempty( AstralMT{jAster} )
                        AstralMT{jAster} = {newMT};
                    else
                        AstralMT{jAster} = { AstralMT{jAster}{:}, newMT};
                    end
                end
            end

            % Store basic objects in object Hierarchy
            % SPBs + Astral MTs stored in an AsterMT
            % SpindleMT + 2 Asters stored in a Spindle
            for jAster = 1 : 2
                if ~isempty( AstralMT{jAster} )
                    AsterObjects{jAster} = AsterMT( dim, imageIn, SPB{jAster}, AstralMT{jAster}{:} );
                else
                    AsterObjects{jAster} = AsterMT( dim, imageIn, SPB{jAster} );
                end

            end
            spindleObj = Spindle( dim, imageIn, {spindleMT, AsterObjects{:} }, props.Spindle);
            spindleObj.findEnvironmentalConditions();
            spindleObj.syncFeaturesWithMap();
            spindleObj.fillParams();
%             spindleObj.displayFeature();
            
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
            spindleMinIntensity = 0.6;
            linewidth = 3;
            brightestPixelAsSPB = 0;

            % Useful variables  
            % zAnisotropy = ceil( sizeVoxelsZ / sizeVoxelsX );
            numVoxelsX = size( imageIn, 1); numVoxelsY = size( imageIn, 2); numVoxelsZ = size( imageIn, 3);

            % Find Strong Signal Regions{{{
            % Convolve the image with a 3D gaussian to bring out the signal from the SPB
            image3DConv = imgaussfilt3( imageIn, 0.5, 'FilterDomain', 'spatial');
            imPlane = mat2gray( max( image3DConv, [], 3) );

            % Keep the strongest signal pixels
            threshOtsu = thresholdOtsu( imPlane( imPlane > 0) );
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
            % dispImg( imPlane); hold on; title('How linescan is done with a line width'); set(gca, 'FontSize', 16)
            % plot( [coords(2,1) coords(2,end)], [coords(1,1) coords(1,end)], 'r-', 'LineWidth', 1)
            for jPix = 1 : length(IntSpindle)
                coordsPerp = round( coords(:, jPix) +  perpMatrix);
                coordsPerp( 1, find( coordsPerp(1,:) > numVoxelsY) ) = numVoxelsY; coordsPerp( 1, find( coordsPerp(1,:) < 1) ) = 1;
                coordsPerp( 2, find( coordsPerp(2,:) > numVoxelsX) ) = numVoxelsX; coordsPerp( 2, find( coordsPerp(2,:) < 1) ) = 1;
            %     plot( coordsPerp(2,:), coordsPerp(1,:), 'c-', 'LineWidth', 1)
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
        function kcBank = findFeaturesDeNovo_KC( image2Find, imageRef)
            % Mitotic Cell: Find Kinetochores
            
            plotFlag_estimate = 0;

            image2Find_Org = image2Find;
            if nargin == 2
%                 image2Find = image2Find .* mat2gray(imageRef);
            end

            [ image2D, idxMax] = max( image2Find, [], 3);
%             image2D = imgaussfilt( image2D, 1):

            % How to find Kinetochores?
            % Iterative thresholding
            threshVal = median( image2D(:) );
            unmaskReg = sum( image2D(:) > threshVal);
            unmaskThresh = 25; % threshold until 35 pixels remain
            while unmaskReg > unmaskThresh
                unmaskReg = sum( image2D(:) > threshVal);
                threshVal = threshVal+2;
            end
            imgThreshed = image2D; imgThreshed( imgThreshed < threshVal) = 0;
            imgThreshed = imgThreshed .* bwconvhull( logical(imgThreshed) );
            
            % Find maxima in the unmasked regions
            maxima = imregionalmax( imgThreshed);
            [y, x] = find( maxima );
            
            % for each (x,y) pair find the z-pixel that has the max intensity, and find the amplitude
            z = 0*x; amplitude= 0*x;
            z = idxMax( sub2ind( size(image2Find), y, x) );
            amplitude = image2Find_Org( sub2ind( size(image2Find), round(y), round(x), round(z) ) );

            % Create the Kinetochores
            disp(sprintf( 'number of kinetochores found = %d', length(x) ) )
            if length(x) > 6
                error('too many kinetochores found. maybe you should fix the spot finder')
            end
            
            dim = length( size( image2Find) ); 
            if dim==2, sigma=[1.2 1.2]; elseif dim==3, sigma=[1.2 1.2 1.0]; end
            for jspot = 1 : length(x)
                kcBank{ jspot} = Kinetochore( [ x(jspot), y(jspot), z(jspot)], amplitude(jspot), sigma, dim, image2Find);
            end
            % Create a Kinetochore Bank for handling and storage
            kcBank = KinetochoreBank( kcBank{:} );

            if plotFlag_estimate
                f = figure;
                ax = axes;
                imagesc( max( image2Find, [], 3) ); colormap gray; axis equal; hold on
                plot( x, y, 'Color', 'g', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2);
            end

        end
        % }}}

        % getUpperLowerBounds {{{
        function [ub, lb] = getUpperLowerBounds( vec, vecLabels, image)
            % Get upper and lower bounds for fitting
            
            ub = vec; lb = vec;
            minVox = 1;
            maxVox = size( image);
            estBkg = median( image( image> 0) );
            maxBkg = max( image(:) );
            minBkg = min( image(:) );
            maxAmp = max( image(:) );
            minAmp = 0; % min SnR is half of max SnR of image
            minSig = [1.2, 1.2, 1.0];
            maxSig = [1.5, 1.5, 1.2];
            

            % find the correct label in vecLabels, and start to place in the correct bounds in the correct places

            % Find index of parameters 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'amplitude') ) );
            idxSig = find( ~cellfun( @isempty, strfind( vecLabels, 'sigma') ) );
            idxP0 = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            idxP1 = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            idxP = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            idxBkg = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );

            % Store upper and lower bounds correctly
            if ~isempty( idxAmp), 
                ub( idxAmp) = maxAmp;
                lb( idxAmp) = minAmp; end
            if ~isempty( idxSig), 
                nF = length( idxSig)/3;
                ub( idxSig) = repmat( maxSig, 1, nF);
                lb( idxSig) = repmat( minSig, 1, nF); end
            if ~isempty( idxP0), 
                nF = length( idxP0)/3;
                ub( idxP0 ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP0 ) = minVox; end
            if ~isempty( idxP1), 
                nF = length( idxP1)/3;
                ub( idxP1 ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP1 ) = minVox; end
            if ~isempty( idxP), 
                nF = length( idxP)/3;
                ub( idxP ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP ) = minVox; end
            if ~isempty( idxBkg), 
                ub( idxP ) = maxBkg;
                lb( idxP ) = minBkg; end

            if any( lb > ub)
                disp( ub )
                disp( lb )
                error('bounds are wrong')
            end

        end
        % }}}

        % getExplorationSpeedVector {{{
        function speedVec = getExplorationSpeedVector( vecLabels)
            % Creates a weighing vector to allow a user to assign different importance to different kinds of parameters. This will infact allow the fitting optimization engine to explore the parameters space at different speeds
            
            % Exlporation Speed : unassigned speeds are kept at 1.0
%             speedAmp = 100;
%             speedBkg = 10;
%             speedSigma = 1;
%             speedPos = 10;
            speedAmp = 1000;
            speedBkg = 10;
            speedSigma = 1;
            speedPos = 1;
            speedVec = ones( size(vecLabels) );

            % Find the index of these speeds
            % Amplitude
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'amplitude') ) );
            speedVec( idxAmp) = speedAmp;

            % Background 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );
            speedVec( idxAmp) = speedBkg;

            % Sigma 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'sigma') ) );
            speedVec( idxAmp) = speedSigma;
            
            % Position
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            speedVec( idxAmp) = speedPos;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            speedVec( idxAmp) = speedPos;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            speedVec( idxAmp) = speedPos;

        end
        % }}}

    end

end
