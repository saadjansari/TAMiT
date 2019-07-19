classdef AsterMT < Organizer
% This is a Microtubule Aster. This sits inside the spindle. It contains 1 SPB (Spindle Pole Body) and a number of microtubules
    properties
    end

    methods
       
        % AsterMT {{{
        function obj = AsterMT( dim, image, spb, varargin)

            if nargin < 3 
                error( 'MT_Array: number of inputs is inconsistent')
            end

            obj = obj@Organizer( dim, image, { spb, varargin{:} }, 'AsterMT');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, propsSPB, propsMT)

            if nargin < 3
                propsMT = {'endPosition', 'amplitude', 'sigma'};
            end
            if nargin < 2
                propsSPB = {'position', 'amplitude', 'sigma'};
            end
            
            vec = [];
            vecLabels = {};

            % get vectors from features

            % get SPB vector. Keep all its values
            [vec_spb, vecLabels_spb] = getVec( obj.featureList{ 1}, propsSPB);
            vec = [vec , vec_spb];
            vecLabels_spb = strcat( 'SPB_', vecLabels_spb);
            vecLabels = { vecLabels{:}, vecLabels_spb{:} };

            % Loop over microtubules and get their vectors. We will truncate its vector to remove the start position
            for jmt = 2 : obj.numFeatures
                [vec_mt, vecLabels_mt] = getVec( obj.featureList{ jmt}, propsMT );
                idxRm = find( strcmp( vecLabels_mt, 'startPosition') );
                vec_mt( idxRm) = [];
                vecLabels_mt( idxRm) = [];
                vec = [vec , vec_mt];
                vecLabels_mt = strcat( 'MT', num2str( jmt-1), '_', vecLabels_mt);
                vecLabels = { vecLabels{:}, vecLabels_mt{:} };
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)

            % Take the vector and find the indexes associated with the SPB parameters
            idxSPB = find( ~cellfun( @isempty, strfind( vecLabels, 'SPB_') ) );
            idxSPBPosition = find( ~cellfun( @isempty, strfind( vecLabels, 'SPB_position') ) );
            % Get the vector for the spindle by removing the spindle substring. Absorb the vector
            vecSPB = vec( idxSPB);
            vecLabelsSPB = erase( vecLabels( idxSPB), 'SPB_');
            obj.featureList{ 1} = absorbVec( obj.featureList{ 1}, vecSPB, vecLabelsSPB );

            % Get the vector for each microtubule, and ask the microtubule objects to absorb the vector
            for jmt = 2 : obj.numFeatures
                
                strMT = ['MT', num2str(jmt-1), '_'];
                idxMT = find( ~cellfun( @isempty, strfind( vecLabels, strMT ) ) );
                vecMT = vec( idxMT );
                vecLabelsMT = erase( vecLabels( idxMT), strMT );
                % Append the parameters for the start position for this microtubules
                vecMT = [ vec( idxSPBPosition) , vecMT ];
                vecLabelsMT_start = repmat( {'startPosition'}, 1, length( idxSPBPosition) );
                vecLabelsMT = { vecLabelsMT_start{:}, vecLabelsMT{:} };
                obj.featureList{ jmt} = absorbVec( obj.featureList{ jmt}, vecMT, vecLabelsMT);

            end
        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures( obj )

            % Get the SPB  position and update the MTs  associated with the SPB
            for jmt = 1 : obj.numMicrotubules
                obj.featureList{1+jmt}.startPosition = obj.featureList{1}.position;
            end

        end
        % }}}

        % findBestMissingFeature {{{
        function missingFeature = findBestMissingFeature( obj, Image2Find, xAngle, xRange)

            % Find Possible missing features
            missingFeatures = findMissingFeatures( obj, Image2Find, xAngle, xRange);

            % Find the best missing feature
            residual = [];
            for jmt = 1 : length( missingFeatures)
                % Find the best MT ( find residual under each MT and pick the biggest one)
                lineAmp{jmt} = Cell.findAmplitudeAlongLine( Image2Find, missingFeatures{jmt}.startPosition, missingFeatures{jmt}.endPosition);
                residual(jmt) = sum( lineAmp{ jmt} );
            end

            % best feature
            [~,idxKeep] = min( residual);
            missingFeature = missingFeatures{idxKeep};
            missingFeature.amplitude = mean( lineAmp{idxKeep} );
            missingFeature.residual = residual(idxKeep);

        end
        % }}}

        % findMissingFeatures {{{
        function missingFeatures = findMissingFeatures( obj, Image2Find, xAngle, xRange)

            % Minimum length of a feature to consider
            Lmin = 5;

            % Find Astral Microtubules
            missingFeatures = MitoticCell.findAstralMicrotubules( Image2Find, obj.featureList{1}.position, xAngle, xRange);

            % Apply the length cutoff
            L = cellfun( @(x) x.length, missingFeatures);
            missingFeatures( find( L<Lmin) ) = [];

            % If missing features aren't found organically, we will guestimate a feature
            if isempty( missingFeatures)
                % Take the image (which is actually the residual from a previous fit), then bleach the SPB region to eliminate really small lines. Find the location of the max residual and draw a feature from the SPB to this point of max residual.
                
                % Bleach the SPB
                Image2Find2D = max( Image2Find, [], 3);
                imgBleached = imgaussfilt( Image2Find2D, 1);
                imgBleached = obj.featureList{1}.BleachSpot( imgBleached, Lmin); 

                [~,indMax] = max( imgBleached(:) );
                [ endY, endX] = ind2sub( size(imgBleached), indMax);
                
                % feature 1
                missingFeatures{1}.startPosition = obj.featureList{1}.position;
                missingFeatures{1}.endPosition = obj.featureList{1}.position;
                missingFeatures{1}.endPosition(1:2) = [ endX, endY];
                lineAmp = Cell.findAmplitudeAlongLine( Image2Find, missingFeatures{1}.startPosition, missingFeatures{1}.endPosition);
                missingFeatures{1}.amplitude = mean( lineAmp);
                missingFeatures{1}.residual = mean( abs( lineAmp - mean(lineAmp) ) );
            end

        end
        % }}}

        % findWorstFeature {{{
        function worstFeature = findWorstFeature( obj, ImageRef)
            % finds the worst microtubule in this aster by looking at the max values under ImageRef.
            % Output :  worstFeature.idx is the idx of the feature (if idx = 0, there are no removable features in this aster)
            %           worstFeature.residual is the summed residual of this feature

            worstFeature.residual = [];
            for jmt = 2 : obj.numFeatures

                % Find the summed residual under each microtubule 
                worstFeature.residual(jmt-1) = sum( Cell.findAmplitudeAlongLine( ImageRef, obj.featureList{jmt}.startPosition, obj.featureList{jmt}.endPosition) );

            end

            % worst microtubule of this aster
            if isempty( worstFeature.residual)
                worstFeature.idx = 0;
            else
                [ worstFeature.residual , worstFeature.idx ] = max( worstFeature.residual);
            end

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.image, [], 3) ); colormap gray; axis equal;
            end

            % Ask subfeatures to display themselves
            % Do this in reverse so that the SPB is plotted last
            for jFeat = obj.numFeatures : -1 : 1
                ax = obj.featureList{jFeat}.displayFeature( ax);
            end

        end
        % }}}

    end

end
