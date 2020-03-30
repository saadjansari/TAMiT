classdef IMTBank < OrganizerMaster
% A IMTBank is a higher level feature that sits inside an interphase cell. It is composed of Microtubule asters 
    properties
        polyOrder % polynomial order of curves
    end

    methods (Access = public)
        
        % IMTBank {{{
        function obj = IMTBank( dim, image, featureList, props2Fit)

            obj = obj@OrganizerMaster( dim, image, featureList, props2Fit, 'IMTBank');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj)

            % Get general environmental properties
            [vec, vecLabels] = getVecEnvironment( obj);

            % Loop over mt arrays and get their vectors. Remove the MTOCs from the fit vectors. 
            for jAster = 1 : obj.numFeatures
                [vec_mtarray, vecLabels_mtarray] = getVec( obj.featureList{ jAster}, obj.props2Fit.fit{obj.dim}.aster.curve);
                vec = [vec , vec_mtarray];
                vecLabels_mtarray = strcat( 'B', num2str( jAster), '_', vecLabels_mtarray);
                vecLabels = { vecLabels{:}, vecLabels_mtarray{:} };
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels);

            % Aster Vector
            % Create the vector for each subfeature MT_Array and send it to the objects for absorption
            for jb = 1 : obj.numFeatures

                strAster = [ 'B', num2str(jb), '_' ];

                % Take the vector and find indices associated with the mt_array. 
                idxA = find( ~cellfun( @isempty, strfind( vecLabels, strAster ) ) );
                vecA = vec( idxA );
                vecLabelsA = erase( vecLabels( idxA ), strAster );
                obj.featureList{ jb} = absorbVec( obj.featureList{ jb}, vecA, vecLabelsA );

            end
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            vecList = {};
            vecLabelsList = {};
            objList = {};

            % bundle vectors
            for jb = 1 : obj.numFeatures
                [ vecList{jb} , vecLabelsList{jb}] = getVec( obj.featureList{jb}, obj.props2Fit.fit{obj.dim}.aster.curve);
                objList{jb} = obj.featureList{jb};
            end
                
            % Prepend environmental parameters to each vec and vecLabel
            [vecE, vecLabelsE] = getVecEnvironment( obj, obj.props2Fit.fit{obj.dim}.Environment);
            
            for jFeat = 1 : length( vecList)
                vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
                vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
            end

        end
        % }}}
        
        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, Image2Find)
            % Searches for missing subfeatures and adds them to the featureList
            % We look for missing microtubule bundles            

            bkg = median( Image2Find( Image2Find(:) > 0) );
            minLength = 20;
            [Image2Find2D, idxZ] = max(Image2Find, [],3);
            if obj.dim == 2
                sigma=[1.2 1.2];
            elseif obj.dim == 3
                sigma=[1.2 1.2 1.0];
            end

            % Run Steerable filter to detect curves
            [imGauss, imSteer] = Methods.FilterImageForCurveDetection( Image2Find2D);
            
            % Remove intensity for curves
            imCurr = obj.simulateFeature( size(Image2Find) );
            if obj.dim == 3
                imCurrMask = max(imCurr,[],3) < 0.1*max(imCurr(:));
            elseif obj.dim == 2
                imCurrMask = imCurr < 0.1*max(imCurr(:));
            end            
            % Get image stats
            stats = Methods.GetImageStats( obj.image, 0);
            
            % Find a curve in the 2D steerable image
            [coord, ~, out] = Methods.FindCurve2D( imSteer.*imCurrMask, 'ThreshInt', stats.ThreshHigh);
            
            % Return if useful curve not found
            if ~out.success || out.Length < minLength 
                successAdd = 0;
                return
            end

            % Get coefficients
            coeffs = Bundle.estimatePolyCoefficients( coord, obj.polyOrder(1:2));
            % coeff Z
            [cX,cY] = Methods.InterpolateCoords( coord(1,:), coord(2,:), 20); 
            indXY = sub2ind( size(Image2Find2D), round(cY), round(cX) );
            [IntSeed, IdxSeed] = max( Image2Find2D( indXY) );
            [ Seed(2), Seed(1)] = ind2sub( size(Image2Find2D), indXY(IdxSeed) );
            % Return if curve seed is not bright
            if max( obj.image( Seed(2), Seed(1), :) ) < stats.ThreshLow
                successAdd = 0;
                return
            end
            % add 3D coefficient 
            if obj.dim == 3
                coeffs{3}(1:2) = [ 0 idxZ( Seed(2), Seed(1) )];
                    coeffs{3}(1) = -0.5;
                if coeffs{3}(2) == 1
                    coeffs{3}(1) = 0.5;
                elseif coeffs{3}(2) == size(Image2Find, 3) 
                    coeffs{3}(1) = -0.5;
                end
            end

            % Find overlap region by thresholding amplitudes
            % Get amplitude along bundle
            amps = Cell.findAmplitudeAlongCurve( max(Image2Find2D,[],3), coeffs ) - bkg;
            amps = smooth( amps, 3);
            amps( amps < bkg) = bkg;
            % Threshold amplitude
            thr = multithresh( amps(:), 1);
            
            % Find parameters vals when amplitude is above threshold
            idx1 = find(amps > thr, 1, 'first');
            idx2 = find(amps > thr, 1, 'last');
            if idx1 == 0
                disp('overlap region is at the start of the curve')
            end
            if idx2 == length(amps) 
                disp('overlap region is ar the end of the curve')
            end
            t = [idx1 idx2]/length(amps);
            t(1) = max( [ t(1) 0.03]);
            t(1) = min( [ t(1) 0.97]);
            t(2) = max( [ t(2) 0.03]);
            t(2) = min( [ t(2) 0.97]);
            
            if abs(diff(t)) < 0.03
                t(1) = max( [ 0.03, t(1)-0.02]);
                t(2) = min( [ 0.97, t(1)+0.02]);
                disp('modified t because values were too close')
            end
            
            % Find amplitude and amplitude enhancement factor
            amp = median( amps(amps <= thr) );
            ef = median( amps(amps > thr) )/amp;
            if ef < 1.5
                disp('enhancement factor is less than 1. set to 1.6')
                ef = 1.6;
            end
            if ef > 4
                disp('enhancement factor is greater than 4. set to 3.9')
                ef = 3.9;
            end
                
            % Create 
            bundle = Bundle( coeffs, amp, sigma, obj.dim, obj.props2Fit.fit{obj.dim}.curve, obj.props2Fit.graphics.curve, t, ef);
            
            if bundle.GetLength() <= minLength
                successAdd = 0;
                return
            end

            % Add the bundle 
            obj.addFeatureToList( bundle);
            % Add feature ID and update the featureMap
            bundle.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            successAdd = 1;
                            
        end
        % }}}
        
        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % For an interphase array, we can only remove  interphase microtubules
            % For each aster we will find a possible microtubule to remove, then we will pick the best of the possibilities and remove it from our feature list.
            % Look at the residual under each astral microtubule
            % Remove the one with the biggest unit residual per length.
            % If there are no microtubules to remove
            successRemove = 1;        

            % If no features to remove, exit the function
            if obj.numFeatures == 0
                successRemove = 0; return
            end

            % Find prominence of features
            p = zeros(1, obj.numFeatures);
            for idx = 1 : obj.numFeatures
                
                % find Prominence
                p(idx) = obj.featureList{ idx }.amplitude * ...
                    obj.featureList{ idx }.GetLength() * ...
                    mean( obj.featureList{ idx }.sigma);

            end
            
            % Find least prominent worst feature
            [~, idxRm ] = min(p);

            % Remove the worst microtubule
            obj.removeFeatureFromList( idxRm);
                
        end
        % }}}
        
        % Update3DFrom2DSpecific {{{
        function obj = Update3DFrom2DSpecific(obj, obj2D)
            % Update background and numAsters

            obj.findEnvironmentalConditions();
            
            % For each aster, set the z value of the seed and the
            % z coefficient of the curves
            for jf = 1 : obj.numFeatures
                seed = obj.featureList{jF}.featureList{1};
                [~,seed.position(3)] = max( obj.image( seed.position(2), seed.position(1), :) );
                seed.sigma(3) = 1.0;
                for jc = 2 : obj.featureList{jf}.numFeatures
                    obj.featureList{jf}.featureList{jc}.startPosition(3) = seed.position(3);
                    obj.featureList{jf}.featureList{jc}.sigma(3) = 1.0;
                    obj.featureList{jf}.featureList{jc}.cZ = 0.0;
                end
            end
                
        end
        % }}}
        
%         function obj = finalizeAddedFeatures(obj)
%             
%             % Reset t parameter for bundles
%             for jb = 1 : obj.numFeatures
%                 obj.featureList{jb}.resetT();
%             end
%             
%         end
        
        % GetStructBasicFeatures {{{
        function feats = GetStructBasicFeatures(obj)

            for jf = 1 : obj.numFeatures
                next = 1;
                feats(jf).type = obj.featureList{jf}.type;
                switch feats(jf).type
                    case 'Spot'
                        feats(next).position = obj.featureList{jf}.position;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        next = next+1;

                    case 'Line'
                        feats(next).startPosition = obj.featureList{jf}.startPosition;
                        feats(next).endPosition = obj.featureList{jf}.endPosition;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        feats(next).length = obj.featureList{jf}.GetLength();
                        feats(next).orientation = obj.featureList{jf}.GetOrientation();
                        next = next+1;

                    case 'Curve'
                        feats(next).startPosition = obj.featureList{jf}.startPosition;
                        feats(next).amplitude = obj.featureList{jf}.amplitude;
                        feats(next).sigma = obj.featureList{jf}.sigma;
                        feats(next).length = obj.featureList{jf}.GetLength();
                        feats(next).orientation = obj.featureList{jf}.GetOrientation();
                        coords = obj.featureList{jf}.GetCoords(); 
                        feats(jf).endPosition = coords(:,end);
                        next = next+1;

                    case 'AsterMT'
                end
            end
        end
        % }}}

    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'IMTBank')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            for jFeat = 1 : length( S.featureList)
                featureList{ jFeat} = AsterMT.loadFromStruct( S.featureList{ jFeat} ); 
            end

            obj = IMTBank( S.dim, S.image, featureList, S.props2Fit);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap();

        end
        % }}}
        
        % findIntBank{{{
        function intBankObj = findIntBank( imageIn, params, props)
            % Finds an interphase MT bank in the given image

            % Get dimensionality and background
            dim = length( size(imageIn) );
            bkg = median( imageIn( imageIn(:) > 0) );

            % Define sigma for this bank
            if dim == 2
                sigma=[1.2 1.2];
            elseif dim == 3
                sigma=[1.2 1.2 1.0];
            end

            % Filter image 
            [imG, ~] = Methods.FilterImageForCurveDetection( imageIn);

            % Find the curves
            vars = Methods.struct2cellvars(params);
            coords = Methods.FindCurves( imG, 'Plot', 0, vars{:}); 
            nBundles = length( coords);
            
            % Create interphase bundles
            for jb = 1 : nBundles

                % Get coord: unite the two curves
                if length( coords{jb}) == 2
                    coord = [ coords{jb}{2}(:, end:-1:2), coords{jb}{1} ];
                else
                    coord = coords{jb}{1};
                end

                % Get coefficients
                coeffs = Bundle.estimatePolyCoefficients( coord, params.PolyOrder);

                % Get amplitude along bundle
                amps = Cell.findAmplitudeAlongCurve( max(imageIn,[],3), coeffs ) - bkg;
                amps = smooth( amps, 3);
                amps( amps < bkg) = bkg;
                % Threshold amplitude
                thr = multithresh( amps(:), 1);
                
                % Find parameters vals when amplitude is above threshold
                idx1 = find(amps > thr, 1, 'first');
                idx2 = find(amps > thr, 1, 'last');
                if idx1 == 0
                    disp('overlap region is at the start of the curve')
                end
                if idx2 == length(amps) 
                    disp('overlap region is ar the end of the curve')
                end
                t = [idx1 idx2]/length(amps);
                t(1) = max( [ t(1) 0.03]);
                t(1) = min( [ t(1) 0.97]);
                t(2) = max( [ t(2) 0.03]);
                t(2) = min( [ t(2) 0.97]);

                % Find amplitude and amplitude enhancement factor
                amp = median( amps(amps <= thr) );
                ef = median( amps(amps > thr) )/amp;
                if ef < 1.5
                    disp('enhancement factor is less than 1. set to 1.6')
                    ef = 1.6;
                end
                if ef > 4
                    disp('enhancement factor is greater than 4. set to 3.9')
                    ef = 3.9;
                end
                    
                % Ensure coefficients are within the image region
                if dim == 3
                    if coeffs{3}(2) == 1
                        coeffs{3}(1) = 0.5;
                    elseif coeffs{3}(2) == size(imageIn, 3)
                        coeffs{3}(1) = -0.5;
                    end
                end
                
                % Create 
                bundles{jb} = Bundle( coeffs, amp, sigma, dim, props.fit{dim}.curve, props.graphics.curve, t, ef);

            end

            % Create the interphase aster bank.
            intBankObj = IMTBank( dim, imageIn, bundles, props);
            intBankObj.polyOrder = params.PolyOrder;
            intBankObj.findEnvironmentalConditions();
            intBankObj.forceInsideMask();

        end
        % }}}

    end
end
