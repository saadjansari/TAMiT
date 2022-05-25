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
        function obj = EstimateFeatures( obj, estimationImage, cTime, cChannel, idxChannel)
        % findFeatures : estimates and finds the features for a monopolar
        % aster
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Estimation
            obj.featureList{ idxChannel, cTime} = obj.EstimateFeaturesNovel( currentFeature, estimationImage);

            % Special Tasks 
            % If cut7 channel, try to get monopolar information from MT channel
            if strcmp( obj.featuresInChannels{ idxChannel}, 'Cut7')
                % Get spindle information from microtubule channel if possible
                % (MICROTUBULE CHANNEL MUST BE ANALYZED FIRST FOR THIS TO WORK)
                mtChannel = find( strcmp( obj.featuresInChannels, 'Microtubule' ) );
                if ~isempty( mtChannel)
                    obj.featureList{ idxChannel, cTime}.harvestMonopolarAsterInfo( obj.featureList{mtChannel, cTime} );
                end
            end
        end
        % }}}

        % EstimateFeaturesNovel {{{
        function feature = EstimateFeaturesNovel( obj, currentFeature, image)

            switch currentFeature
                case 'Microtubule'
                    feature = MonopolarCell.findFeaturesDeNovo_MT( image, obj.params.estimate.monopolar, obj.featureProps.monopolarAster);

                case 'Kinetochore'
                    feature = MonopolarCell.findFeaturesDeNovo_KC( image, obj.params.estimate.kcbank);

                case 'Cut7'
                    feature = MonopolarCell.findFeaturesDeNovo_Cut7( image, obj.params.estimate.cut7dist);
            end

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function MonopolarObj = findFeaturesDeNovo_MT( imageIn, params, props)
            % Monopolar Cell:
            %   Find the Pole
            %   From astral microtubules

            % #################################
            % ###########  Params  ############
            % #################################
            
            displayFlag = params.display;
            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            sigma = params.common.sigma(1:dim);
            bkg = median( imageIn( imageIn(:) > 0) );
            astralMTRepresentation = params.astralMTRepresentation;
            astralMinLength = params.astralMinLength;

            % #################################
            % #########  Spindle Pole #########
            % #################################

            % Set up parameters for pole finder
            params_pole_finder.poleMinIntensity = params.poleMinIntensity;
            params_pole_finder.poleDeterminationSensitivity = params.poleDeterminationSensitivity;
            params_pole_finder.visuals = params.visuals;
            params_pole_finder.visuals_path = params.visuals_path;
            params_pole_finder.verbose = params.common.verbose;

            % Extract the pole
            pole = MonopolarCell.findThePole( imageIn, params_pole_finder); 

            % Find amplitude at spindle pole, then subtract background.
            poleAmp = imageIn( round(pole.position(2)), round(pole.position(1)), round(pole.position(3)) ) - bkg;
            
            % Force pole amplitude to be above atleast equal to the
            % background.
            if poleAmp < 0
                poleAmp = bkg;
                if params.common.verbose>1
                    warning(sprintf('\n%s\n%s\n%s\n%s\n',...
                        'POTENTIALLY BAD DETECTION!!!',...
                        'AMPLITUDE OF SPINDLE POLE FOUND WAS LESS THAN BKG!!!',...
                        'FORCING SPINDLE POLE AMPLITUDE TO BE EQUAL TO BKG!!!',...
                        'FINAL RESULTS MAY BE INACCURATE!!!'))
                end
            end
            
            % Create a Spot feature for the Spindle Pole Body
            Pole = Spot( pole.position, poleAmp, sigma, dim, props.fit{dim}.spot, props.graphics.aster.spot); Pole.label = 'spb';
            
            % #################################
            % ##########  Astral MTs  #########
            % #################################

            % Astral MTs
            lines = {};
            
            if params.astralMT

                % Find Astral Microtubules
                Aster.MT = MitoticCell.findAstralMicrotubules( imageIn, pole.position, [], [] );

                for jmt = 1 : length( Aster.MT )
                    
                    % Find median amplitude along line
                    lineAmp = median( Cell.findAmplitudeAlongLine( imageIn, Aster.MT{jmt}.startPosition, Aster.MT{jmt}.endPosition ) )-bkg;
                    
                    % Create Line objects for Astral MTs
                    newMT = Line( Aster.MT{jmt}.startPosition, Aster.MT{jmt}.endPosition, lineAmp, sigma, dim, props.fit{dim}.line, props.graphics.aster.line);
                    newMT.repr = astralMTRepresentation; 
                    newMT.label = 'astralFYMonopolar';

                    % Add lines longer than min allowed length
                    if norm( [Aster.MT{jmt}.startPosition - Aster.MT{jmt}.endPosition]) > astralMinLength
                        lines = { lines{:}, newMT };
                    end
                end
            end

            % ##############################################
            % ##########  Setup Feature Heirarchy  #########
            % ##############################################
            
            % Aster MT Objects
            % SPB + Astral MTs stored in an AsterMT
            AsterObj  = AsterMT( dim, Pole, lines{:} );

            % Monopolar Aster Feature 
            MonopolarObj = MonopolarAster( dim, imageIn, {AsterObj}, props);
            
            % Find enviornment conditions for spindle
            MonopolarObj.findEnvironmentalConditions();
            
        end
        % }}}
        
        % findThePole {{{
        function [pole, h] = findThePole( imageIn, params)
            % Find a bright spindle pole
            % Input Parameters:
            % -----------------
            %   1. imageIn : 3D array (XYZ intensity data)
            %   2. params  : struct
            %       A. poleMinIntensity : float (Minimum intensity of pole
            %           relative to brightest image pixel) (Range 0-1)
            %       B. poleDeterminationSensitivity : float (Extended
            %           Region Sensitivity) (Range 0-1)
            %       C. visuals : boolean (turn on visuals)
            %       D. visuals_path : str (location to save the visuals)
            %       E. verbose : 0,1,2 (print relevant info to log)
            
            imageIn_old = imageIn;
            h=0;
            
            % Extract Params
            % Minimum intensity of pole relative to brightest image pixel
            poleMinIntensity = params.poleMinIntensity;
            
            % Extended Region Sensitivity
            poleDeterminationSensitivity = params.poleDeterminationSensitivity;
            
            % Other parameters
            visuals = params.visuals; % to turn on visuals
            visuals_path = params.visuals_path; % location to save the visuals
            verbose = params.verbose; % descriptions ON
                                    
            % Find the mask
            imMask3D = imageIn > 0;

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

            if verbose>0; fprintf( ' - Pole Mean Intensity = %.2f\n', poleInt); end
            if poleInt < poleMinIntensity
                if verbose>1; warning( 'Mean pole intensity < 70% of the max image intensity'); end
            end

            % Now we can obtain the shape parameters of the spindle
            pole.position = cell2mat( SScell( 1, idx) );
            % }}}

            % find Z coordinates of the pole 
            [~, zIdx] = max( image3DConv, [], 3);
            voxelZ = zIdx( round(pole.position(2)), round(pole.position(1)) );
            pole.position = [ pole.position, voxelZ];

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
