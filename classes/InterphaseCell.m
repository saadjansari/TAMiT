classdef InterphaseCell < Cell
% This is a specialized cell of superclass Cell that is used for an interphase cell.
    properties
    end

    methods
        
        % InterphaseCell {{{
        function obj = InterphaseCell( image, lifetime, species, featuresInChannels, channelsToFit, settings)
        % InterphaseCell: constructor function for InterphaseCell object     

            obj = obj@Cell( image, lifetime, species, featuresInChannels, channelsToFit, 'Interphase', settings);

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

                        obj.featureList{ idxChannel, cTime} = MitoticCell.findFeaturesDeNovo_MT( obj.image(:,:,:,cTime, cChannel), obj.settings.flags.debug );
                    
                    otherwise

                        error('InterphaseCell: unknown feature')

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
        function imtBankObj = findFeaturesDeNovo_MT( imageIn, displayFlag )

            % Interphase Cell:
            %   Find iMTOCs (interphase microtubule organizing centers)
            %   From each iMTOC, find micrtoubule curves 

            displayFlag = displayFlag;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);

            props.MTOC= {'position', 'amplitude', 'sigma'};
            props.AsterMT = {'endPolyCoef', 'amplitude', 'sigma'};
            props.iMT = {'polyCoef', 'amplitude', 'sigma'};
            props.iMTBank = {'background'};
            display.MTOC= {'Color', [0.7 0 0.7] , 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2};
            display.iMT= {'Color', [1 0 1] , 'LineWidth', 2};
            display.AsterMT = {'Color', [1 0 1] , 'LineWidth', 2};
            if dim==3, sigma=[1.2 1.2 1.0]; elseif dim==2, sigma=[1.2 1.2]; end
            bkg = median( imageIn( imageIn(:) > 0) );

            % Find the iMTOCs 
            
            % Find the Spindle. 
            if displayFlag
                f = figure;
                ax = axes;
                imagesc( max(imageIn, [], 3) ); colormap gray; axis equal; hold on
                [mtocs, ax]  = MitoticCell.findIMTOCs( imageIn, ax);
            else
                mtocs = InterphaseCell.findIMTOCs( imageIn);
            end

            % Create the MTOC objects
            for jMTOC = 1 : mtocs.num
                
                % get amplitude of mtoc
                mtoc_amp = mtocs.amp(jMTOC);
                if mtoc_amp < 0
                    mtocs_amp = bkg;
                    warning( 'findFeaturesMT_deNovo : forcing MTOCamp to be > 0')
                end

                MTOC{ jMTOC} = Spot( mtocs.position, mtoc_amp, sigma, dim, props.MTOC, display.MTOC);
            end

            % Find interphase MTs
            for jMTOC = 1 : length( MTOC)
                
                if displayFlag 
                    [imt{jMTOC}, ax] = InterphaseCell.findIMTs( imageIn, MTOC{ jMTOC}.position, ax );
                else
                    imt{jMTOC} = InterphaseCell.findIMTs( imageIn, MTOC{ jMTOC}.position );
                end

                % Create interphase MTs
                iMT = {};
                for jmt = 1 : length( iMT{jMTOC} )
                    
                    curveAmp = median( Cell.findAmplitudeAlongCurve( imageIn, iMT{jMTOC}{jmt}.polyCoef ) ) - bkg;
                    iMT{jmt} = Curve( MTOC{jMTOC}.position, iMT{jMTOC}{jmt}.polyCoef, curveAmp, sigma, dim, props.iMT, display.iMT);
                    iMT{jmt}.findVoxelsInsideMask( logical(imageIn) );

                end

                % Store iMTOC + iMTs in an AsterMT
                AsterObjects{ jMTOC} = AsterMT( dim, MTOC{ jMTOC}, iMT{:} );

            end

            % AsterMTs stored in a iMTBank 
            iMTBankObj = IMTBank( dim, imageIn, AsterObjects, props.iMTBank);
            iMTBankObj.findEnvironmentalConditions();
            iMTBankObj.syncFeaturesWithMap();
            
        end
        % }}}
        % findIMTOCs {{{
        function [mtocs, ax] = findIMTOCs( imageIn, ax)
            
            if nargin < 2 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end

            % Params
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
        function [ imts, ax] = findIMTs( imageIn, startPoint, ax)
            % Finds 2 opposite curved lines directed away from the start point

            if nargin < 3 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end
            
            Image2Find = imgaussfilt( imageIn, 1);

            if nargout > 1 
                [iMTs, ax] = Cell.find_manyCurves( Image2Find, startPoint, ax);
            else
                iMTs = Cell.find_manyCurves( Image2Find, startPoint);
            end

            if isempty( iMTs )
                nomt = 1;
            else
                nomt = 0; 
            end

        end
        % }}}
        
        % }}}

    end

end
