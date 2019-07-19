classdef MT_Array < Feature
% This is a Microtubule Aster. This sits inside the spindle. It contains 1 SPB (Spindle Pole Body) and a number of microtubules
    properties
        numMicrotubules 
        featureList
    end

    methods
       
        % MT_Array {{{
        function obj = MT_Array( spb, varargin)

           if nargin == 0 
                error( 'MT_Array: number of inputs is inconsistent')
            end

            obj = obj@Feature( spb.position, spb.amplitude, spb.sigma, spb.dim, spb.ref_image, 'MT_Array');

            if nargin == 1
                obj.featureList = { spb };
            elseif nargin >= 2
                obj.featureList = { spb, varargin{:} };
            end

            obj.numMicrotubules = length( varargin );

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props)

            vec = [];
            vecLabels = {};

            % get vectors from features

            % get SPB vector. Keep all its values
            [vec_spb, vecLabels_spb] = getVec( obj.featureList{ 1}, props.spb );
            vec = [vec , vec_spb];
            vecLabels_spb = strcat( 'SPB_', vecLabels_spb);
            vecLabels = { vecLabels{:}, vecLabels_spb{:} };

            % Loop over microtubules and get their vectors. We will truncate its vector to remove the start position
            for jmt = 1 : length( obj.featureList)-1
                [vec_mt, vecLabels_mt] = getVec( obj.featureList{ 1+jmt}, props.mt );
                idxRm = find( strcmp( vecLabels_mt, 'startPosition') );
                vec_mt( idxRm) = [];
                vecLabels_mt( idxRm) = [];
                vec = [vec , vec_mt];
                vecLabels_mt = strcat( 'MT', num2str( jmt), '_', vecLabels_mt);
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
            for jmt = 1 : obj.numMicrotubules
                
                strMT = ['MT', num2str(jmt), '_'];
                idxMT = find( ~cellfun( @isempty, strfind( vecLabels, strMT ) ) );
                vecMT = vec( idxMT );
                vecLabelsMT = erase( vecLabels( idxMT), strMT );
                % Append the parameters for the start position for this microtubules
                vecMT = [ vec( idxSPBPosition) , vecMT ];
                vecLabelsMT_start = repmat( {'startPosition'}, 1, length( idxSPBPosition) );
                vecLabelsMT = { vecLabelsMT_start{:}, vecLabelsMT{:} };
                obj.featureList{ 1+jmt} = absorbVec( obj.featureList{ 1+jmt}, vecMT, vecLabelsMT);

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

        % simulateFeature {{{ 
        function imageOut = simulateFeature( obj, imageIn)

            if nargin < 2
                imageIn = obj.refImage;
            end

            % Draw the Spindle pole Body
            imageOut = simulateFeature( obj.featureList{1}, imageIn);

            % Draw the microtubules
            for jmt = 1 : obj.numMicrotubules
                imageOut = simulateFeature( obj.featureList{1+jmt}, imageOut);
            end

        end
        % }}}

        % addFeature {{{
        function obj = addFeature( obj, MT)

            % Check to ensure MT is an object of type microtubule
            if ~strcmp( MT.type, 'MT_linear')
                error('addFeature: feature object to add must be of the correct type'); end

            % Add MT
            obj.featureList = { obj.featureList{:}, MT };
            obj.numMicrotubules = obj.numMicrotubules + 1;

        end
        % }}}

        % removeFeature {{{
        function obj = removeFeature( obj, idxFeature) 

            obj.featureList( idxFeature) = [];
            obj.numMicrotubules = obj.numMicrotubules - 1;

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
                imgBleached = obj.featureList{1}.BleachIt( imgBleached, Lmin); 

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
            for jmt = 1 : obj.numMicrotubules

                % Find the summed residual under each microtubule 
                worstFeature.residual(jmt) = sum( Cell.findAmplitudeAlongLine( ImageRef, obj.featureList{1+jmt}.startPosition, obj.featureList{1+jmt}.endPosition) );

            end

            % worst microtubule of this aster
            if isempty( worstFeature.residual)
                worstFeature.idx = 0;
            else
                [ worstFeature.residual , worstFeature.idx ] = max( worstFeature.residual);
            end

        end
        % }}}

        % makeCopyDeep {{{
        function objCopy = makeCopyDeep( obj)
            % makes a deep copy of the handle object

            % first make a shallow copy
            objCopy = copy( obj);

            % now make it deep by copying obj.featureList
            for jFeat = 1 : length( obj.featureList)
                objCopy.featureList{jFeat} = copy( obj.featureList{jFeat} ); % copies the subfeatures
            end

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.ref_image, [], 3) ); colormap gray; axis equal;
            end

            % Ask subfeatures to display themselves
            for jFeat = 1 : 1 + obj.numMicrotubules
                ax = obj.featureList{jFeat}.displayFeature( ax);
            end

        end
        % }}}

    end
end
