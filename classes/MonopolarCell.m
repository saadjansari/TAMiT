classdef MonopolarCell < Cell
% This is a specialized cell of superclass Cell that is used for a monopolar spindle cell.
    properties
        phase = 'monopolar'% cell phase (e.g. prophase, metaphase, anaphase) if available
    end

    methods
        
        % MonopolarCell {{{
        function obj = MonopolarCell( imageData, featuresInChannels, channelsToFit, params, varargin)
        % MitoticCell : constructor function for MitoticCell object     

            obj = obj@Cell( imageData, featuresInChannels, channelsToFit, params, varargin{:}, 'Type', 'Monopolar');

        end
        % }}}

        % EstimateFeatures {{{
        function obj = EstimateFeatures( obj, cTime, cChannel, idxChannel)
        % findFeatures : estimates and finds the features 
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Get Image for estimation
            estimationImage = obj.imageData.GetImage;
            estimationImage = estimationImage( :,:,:,cTime, cChannel);

            % Get Start time 
            lifetime = obj.imageData.GetLifetime;
            startTime = lifetime(1);

            % Novel Estimation for first frame
            if cTime == startTime
                obj.featureList{ idxChannel, cTime} = obj.EstimateFeaturesNovel( currentFeature, estimationImage);
            % Propagate old feature for later frames
            else 
                obj = obj.PropagateOldFeature( idxChannel, cTime);
            end

            % Special Tasks 
            % If cut7 channel, try to get spindle information from MT channel
            %if strcmp( obj.featuresInChannels{ idxChannel}, 'Cut7')
                %% Get spindle information from microtubule channel if possible
                %% (MICROTUBULE CHANNEL MUST BE ANALYZED FIRST FOR THIS TO WORK)
                %mtChannel = find( strcmp( obj.featuresInChannels, 'Microtubule' ) );
                %if ~isempty( mtChannel)
                    %obj.featureList{ idxChannel, cTime}.harvestSpindleInfo( obj.featureList{mtChannel, cTime} );
                %end

            %end

        end
        % }}}

        % EstimateFeaturesNovel {{{
        function feature = EstimateFeaturesNovel( obj, currentFeature, image)
            disp('- DeNovo') 

            switch currentFeature
                case 'Microtubule'
                    feature = MonopolarCell.findFeaturesDeNovo_MT( image, obj.params.estimate.monopolar);

                case 'Kinetochore'
                    feature = MonopolarCell.findFeaturesDeNovo_KC( image, obj.params.estimate.kcbank);

                case 'Cut7'
                    feature = MonopolarCell.findFeaturesDeNovo_Cut7( image, obj.params.estimate.cut7dist);
            end

        end
        % }}}

        % PropagateOldFeature {{{
        function obj = PropagateOldFeature(obj, cChannel, cTime)

            disp('- Propagating old feature') 

            % Find the most recent good frame
            bestFrame = cTime-1;
            while obj.featureList{ cChannel, bestFrame} ~= obj.featureList{ cChannel, bestFrame}
                bestFrame = bestFrame - 1;
            end

            % Duplicate feature from best recent frame
            obj.featureList{ cChannel, cTime} = obj.featureList{ cChannel, bestFrame}.copyDeep();

            % Ensure the image is from the actual frame
            Image = obj.imageData.GetImage();
            obj.featureList{ cChannel, cTime}.image = Image(:,:,:, cTime, cChannel);

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function MonopolarObj = findFeaturesDeNovo_MT( imageIn, params)

            % Monopolar Cell:
            %   Find astral microtubules from single pole

            displayFlag = params.display;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);

            props.Pole= {'position', 'amplitude', 'sigma'};
            props.AsterMT = {'endPosition', 'amplitude', 'sigma'};
            props.Aster = {'background'};
            display.Pole = {'Color', [0.7 0 0.7] , 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2};
            display.AsterMT = {'Color', [1 0 1] , 'LineWidth', 3};
            if dim==3, sigma=[1.2 1.2 1.0]; elseif dim==2, sigma=[1.2 1.2]; end
            bkg = median( imageIn( imageIn(:) > 0) );

            % Find the Pole
            pole = MonopolarCell.findThePole( imageIn); 

            % Create Pole
            if params.pole

                % Find amplitude
                poleAmp = imageIn( round(pole.position(2)), round(pole.position(1)), round(pole.position(3)) ) - bkg;
                if poleAmp < 0;
                    poleAmp = bkg;
                    warning( 'findFeaturesMT_deNovo : forcing PoleAmp to be > 0')
                end
                Pole = Spot( pole.position, poleAmp, sigma, dim, props.Pole, display.Pole);
            else
                error('Monopolar Cell must have a microtubule pole')
            end

            % Astral MTs
            lines = {};
            if params.mt

                % Find Astral Microtubules
                Aster.MT = MitoticCell.findAstralMicrotubules( imageIn, pole.position, [], [] );

                % Create Astral MTs
                for jmt = 1 : length( Aster.MT )

                    lineAmp = median( Cell.findAmplitudeAlongLine( imageIn, Aster.MT{jmt}.startPosition, Aster.MT{jmt}.endPosition ) )-bkg;

                    newMT = Line( Aster.MT{jmt}.startPosition, Aster.MT{jmt}.endPosition, lineAmp, sigma, dim, props.AsterMT, display.AsterMT);

                    newMT.findVoxelsInsideMask( logical(imageIn) );
                    lines = { lines{:}, newMT };

                end

            end

            % Store basic objects in object Hierarchy
            
            % Aster MT Objects
            % Pole + Astral MTs stored in an AsterMT
            AsterObj  = AsterMT( dim, Pole, lines{:} );

            % Monopolar Aster Feature 
            MonopolarObj = MonopolarAster( dim, imageIn, {AsterObj}, props.Aster);

            if displayFlag
                fName = 'estimate_monopolar';
                f = figure( 'NumberTitle', 'off', 'Name', fName); 
                ax = axes;
                imagesc( ax, max(imageIn, [], 3) ); colormap gray; axis equal; hold on;
                MonopolarObj.displayFeature(ax)
                %sName = [fitInfo.saveDirectory, filesep, sName];
                %export_fig( sName, '-png', '-nocrop', '-a1') 
            end
            MonopolarObj.findEnvironmentalConditions();
            MonopolarObj.syncFeaturesWithMap();
            
        end
        % }}}
        
        % findThePole {{{
        function [pole, ax] = findThePole( imageIn, ax)
            
            if nargin < 2 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end

            % Params
            poleMinIntensity = 0.75;
            poleDeterminationSensitivity = 0.2;
            imMask3D = imageIn > 0;

            % Useful variables  
            numVoxelsX = size( imageIn, 1); numVoxelsY = size( imageIn, 2); numVoxelsZ = size( imageIn, 3);

            % Find Strong Signal Regions{{{
            % Convolve the image with a 3D gaussian to bring out the signal from the SPB
            image3DConv = imgaussfilt3( imageIn, 1, 'FilterDomain', 'spatial') .* imMask3D;
            imPlane = mat2gray( max( image3DConv, [], 3) );

            % Keep the strongest signal pixels
            threshOtsu = max( [thresholdOtsu( imPlane( imPlane > 0) ) 0.7]);
            imPlaneStrong = imPlane;
            imPlaneStrong( imPlaneStrong < threshOtsu) = 0;

            imMask = imextendedmax( imPlaneStrong, poleDeterminationSensitivity);
            imPlaneStrongBool = imPlaneStrong;
            imPlaneStrongBool( imPlaneStrongBool > 0 ) = 1;

            % }}}

            % Pick Brightest Strong Region and find its shape {{{
            % we'll pick the maxima which has the biggest mean intensity
            SS = regionprops( imMask, imPlane, 'Centroid', 'MeanIntensity');
            SScell = struct2cell( SS);
            meanInts = cell2mat( SScell(2,:) );
            [ poleInt, idx] = max( meanInts);

            if poleInt < 0.7
                error( 'The mean pole intensity is less than 70% of the maximum intensity in the image');
            else
                fprintf( 'Success: Pole Mean Intensity = %.2f\n', poleInt) 
            end

            % Now we can obtain the shape parameters of the spindle
            pole.position = cell2mat( SScell( 1, idx) );
            % }}}

            % find Z coordinates of the pole 
            [~, zIdx] = max( image3DConv, [], 3);
            voxelZ = zIdx( round(pole.position(2)), round(pole.position(1)) );
            pole.position = [ pole.position, voxelZ];

            if nargout > 1
                plot( AstersX, AstersY, 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 4, 'Color', 'g' )
            end

        end
        % }}}
        
        % findAstralMicrotubules {{{
        function [AstralMicrotubules, ax] = findAstralMicrotubules( imageIn, startPoint, ax)
            % Finds straight lines directed radially away from the start point
            % Removes the line that matches the spindle angle

            if nargin < 3 && nargout > 1
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
