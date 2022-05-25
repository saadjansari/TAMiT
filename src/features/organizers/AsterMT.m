classdef AsterMT < Organizer
% This is a Microtubule Aster. This sits inside the spindle. It contains 1 SPB (Spindle Pole Body) and a number of microtubules
    properties
    end

    methods
       
        % AsterMT {{{
        function obj = AsterMT( dim, spb, varargin)

            if nargin < 2 
                error( 'MT_Array: number of inputs is inconsistent')
            end

            obj = obj@Organizer( dim, { spb, varargin{:} }, 'AsterMT');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels,ub,lb] = getVec( obj, propsSPB, propsMT)

            if nargin < 3
                propsMT = {'endPosition', 'theta','length','amplitude', 'sigma'};
            end
            if nargin < 2
                propsSPB = {'position', 'amplitude', 'sigma'};
            end
            
            vec = []; lb = []; ub = [];
            vecLabels = {};

            % get vectors from features

            % get SPB vector. Keep all its values
            try
                [vec_spb, vecLabels_spb, ub_spb, lb_spb] = getVec( obj.featureList{ 1}, propsSPB);
            catch
                [vec_spb, vecLabels_spb] = getVec( obj.featureList{ 1}, propsSPB);
            end
            vec = [vec , vec_spb]; lb = [lb, lb_spb]; ub = [ub, ub_spb];
            vecLabels_spb = strcat( 'SPB_', vecLabels_spb);
            vecLabels = { vecLabels{:}, vecLabels_spb{:} };

            % Loop over microtubules and get their vectors. We will truncate its vector to remove the start position
            for jmt = 2 : obj.numFeatures
                try
                    [vec_mt, vecLabels_mt, ub_mt, lb_mt] = getVec( obj.featureList{ jmt}, propsMT );
                catch
                    [vec_mt, vecLabels_mt] = getVec( obj.featureList{ jmt}, propsMT );
                end
                idxRm = find( strcmp( vecLabels_mt, 'startPosition') );
                vec_mt( idxRm) = [];
                vecLabels_mt( idxRm) = [];
                vec = [vec , vec_mt]; lb = [lb, lb_mt]; ub = [ub, ub_mt];
                vecLabels_mt = strcat( 'MT', num2str( jmt-1), '_', vecLabels_mt);
                vecLabels = { vecLabels{:}, vecLabels_mt{:} };
            end
            if isempty(lb)
                clearvars lb ub
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels, errBoolean)
            
            if nargin < 4
                errBoolean = 0;
            end
            
            % Take the vector and find the indexes associated with the SPB parameters
            idxSPB = find( ~cellfun( @isempty, strfind( vecLabels, 'SPB_') ) );
            idxSPBPosition = find( ~cellfun( @isempty, strfind( vecLabels, 'SPB_position') ) );
            % Get the vector for the spindle by removing the spindle substring. Absorb the vector
            vecSPB = vec( idxSPB);
            vecLabelsSPB = erase( vecLabels( idxSPB), 'SPB_');
            obj.featureList{ 1} = absorbVec( obj.featureList{ 1}, vecSPB, vecLabelsSPB, errBoolean );

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
                obj.featureList{ jmt} = absorbVec( obj.featureList{ jmt}, vecMT, vecLabelsMT, errBoolean);

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
        function [missingFeature, succ] = findBestMissingFeature( obj, imOrg, Image2Find, xAngle, xRange)
            % featType = 'Line' or 'Curve'
            
            if nargin <6
                featType = 'Line';
            end
            % Find Possible missing features
            missingFeatures = findMissingFeatures( obj, Image2Find, xAngle, xRange);

            % remove weak lines
            idxRm = []; bkg = median( imOrg( imOrg ~= 0));
            for j1 = 1 : length(missingFeatures)
                if missingFeatures{j1}.amplitude < bkg % i.e SNR is less than 2
                    idxRm = [idxRm, j1];
                end
            end
            missingFeatures( idxRm) = [];
            if isempty(missingFeatures)
                missingFeature = []; succ = 0; return
            end
            
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
            succ = 1;

        end
        % }}}

        % findMissingFeatures {{{
        function missingFeatures = findMissingFeatures( obj, Image2Find, xAngle, xRange)
            % featType = 'Line' or 'Curve'
            
            % Minimum length of a feature to consider
            Lmin = 4;

            % Find Astral Microtubules in residual image
            image2D = max( Image2Find, [], 3);

            % radially integrate in phi
            [ phiIntensity, phiValues] = Cell.radIntegrate2D( image2D, obj.featureList{1}.position,8, 20);
%             figure; plot( phiValues, phiIntensity, 'r-');

            % find peaks in Phi Intensity
            imVals = image2D( image2D ~= 0);
            mtBkg = median( imVals );
            minPkHeight = mtBkg + std( imVals );

            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, 'MinPeakHeight', minPkHeight );
            peakPhi = mod( peakPhi, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )

            % try remove angle belonging to spindle.
            missingFeatures = {};
            if ~isempty(xAngle)
                spindleAngle = mod( xAngle, 2*pi); 
                idxRm = find( abs(peakPhi-spindleAngle) < xRange | abs(peakPhi-spindleAngle+2*pi) < xRange | abs(peakPhi-spindleAngle-2*pi) < xRange);
                peakPhi( idxRm)=[];

                if length(peakPhi) == 0
                    return
                end
            end
            
            % Find position of lines
            for jLine = 1 : length( peakPhi)
                cLine = Cell.find_singleLine3D( imgaussfilt(Image2Find,1), obj.featureList{1}.position, peakPhi( jLine) );
                if ~isempty( cLine)
                    missingFeatures = { missingFeatures{:}, cLine}; 
                end
            end

            % Apply the length cutoff
            L = cellfun( @(x) x.length, missingFeatures);
            missingFeatures( find( L<Lmin) ) = [];
            for j1 = 1 : length(missingFeatures)
                lineAmp = Cell.findAmplitudeAlongLine( max(Image2Find,[],3), missingFeatures{j1}.startPosition(1:2), missingFeatures{j1}.endPosition(1:2));
                missingFeatures{j1}.amplitude = mean( lineAmp( ceil(end/2): end));
                missingFeatures{j1}.residual = mean( abs( lineAmp - mean(lineAmp) ) );
            end
            
            % If missing features aren't found organically, we will guestimate a feature
            if isempty( missingFeatures)
                % Take the image (which is actually the residual from a previous fit), then bleach the SPB region to eliminate really small lines. Find the location of the max residual and draw a feature from the SPB to this point of max residual.
                
                % Bleach the SPB
                Image2Find2D = max( Image2Find, [], 3);
                imgBleached = imgaussfilt( Image2Find2D, 1);
                imgBleached = obj.featureList{1}.BleachSpot( imgBleached, Lmin); 

                conBleach = 1;
                while conBleach
                    
                    [~,indMax] = max( imgBleached(:) );
                    [ endY, endX] = ind2sub( size(imgBleached), indMax);
                    if isempty(xRange) && isempty(xAngle)
                        conBleach = 0;
                    elseif abs( mod( atan2( endY-obj.featureList{1}.position(2), endX-obj.featureList{1}.position(1) ), 2*pi) - xAngle ) > xRange
                        conBleach = 0;
                    else
                        imgBleached( endY-1:endY+1, endX-1:endX+1) = 0;
                    end
                end
                
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

        % findHighResFeature {{{
        function [feat, success] = findHighResFeature( obj, ImageRef)
            % finds the worst microtubule in this aster by looking at the
            % residual( max values under ImageRef)
            % Output :  worstFeature.idx is the idx of the feature (if idx = 0, there are no removable features in this aster)
            %           worstFeature.residual is the summed residual density of this feature

            feat.residual = [];
            feat.idx = 0;
            feat.feature = [];
            success = 1;
            
            % If no features to remove, exit the function
            if obj.numFeatures == 1
                success = 0;
                return
            end
            
            res = zeros(1, obj.numFeatures);
            for jmt = 2 : obj.numFeatures

                % Get length of feature
                len = obj.featureList{jmt}.GetLength();
                
                % Find the summed residual under each microtubule
                switch obj.featureList{jmt}.type
                    case 'Line'
                        res(jmt) = sum( Cell.findAmplitudeAlongLine( ImageRef, obj.featureList{jmt}.startPosition, obj.featureList{jmt}.endPosition) ) / len;
                    case 'Curve'
                        res(jmt) = sum( Cell.findAmplitudeAlongCurve( ImageRef, obj.featureList{jmt}.GetCoeffFull() ) ) /len;
                end
            end

            % worst microtubule of this aster
            [ feat.residual , feat.idx ] = max( res);
            feat.feature = obj.featureList{feat.idx};

        end
        % }}}
        
        % findLeastPromFeature {{{
        function [feat, success] = findLeastPromFeature( obj)
            % finds the worst microtubule in this aster by looking at prominence of features.
            % Prominence = amp * length * mean(sigma)
            
            feat.idx = 0;
            feat.feature = [];
            feat.p = [];
            success = 1;

            % If no features to remove, exit the function
            if obj.numFeatures == 1
                success = 0;
                return
            end

            % Find prominence of features
            p = zeros(1, obj.numFeatures);
            for idx = 2 : obj.numFeatures
                
                % find Prominence
                p(idx) = obj.featureList{ idx }.amplitude * ...
                    obj.featureList{ idx }.GetLength() * ...
                    mean( obj.featureList{ idx }.sigma);

            end
            p(1) = max(p)+1;
            
            % Find least prominent worst feature
            [feat.p, feat.idx ] = min(p);
            feat.feature = obj.featureList{feat.idx};
            
        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, varargin)

            if nargin < 2
                error('must provide an axes to display feature on')
            end

            % Ask subfeatures to display themselves
            % Do this in reverse so that the SPB is plotted last
            for jFeat = obj.numFeatures : -1 : 1
                ax = obj.featureList{jFeat}.displayFeature( ax, varargin{:});
            end

        end
        % }}}

        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;

            for jFeat = 1 : length( obj.featureList)
                S.featureList{ jFeat} = saveAsStruct( obj.featureList{jFeat} );
            end

        end
        % }}}
        
        % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
           % Get 2D projection of feature
            
        end
        % }}}
        

    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'AsterMT')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = Spot.loadFromStruct( S.featureList{ 1} ); 
            for jFeat = 2: length( S.featureList)
                try
                    featureList{ jFeat} = Line.loadFromStruct( S.featureList{jFeat} );
                catch
                    featureList{ jFeat} = Curve.loadFromStruct( S.featureList{jFeat} );
                end
            end

            obj = AsterMT( S.dim, featureList{:});
%             obj = obj@Organizer( S.dim, featureList, S.type);

        end
        % }}}
        
    end

end
