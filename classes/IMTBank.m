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
                [vec_mtarray, vecLabels_mtarray] = getVec( obj.featureList{ jAster}, obj.props2Fit.fit{obj.dim}.aster.spot, ...
                    obj.props2Fit.fit{obj.dim}.aster.curve);
                vec = [vec , vec_mtarray];
                vecLabels_mtarray = strcat( 'A', num2str( jAster), '_', vecLabels_mtarray);
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
            for jAster = 1 : obj.numFeatures

                strAster = [ 'A', num2str(jAster), '_' ];

                % Take the vector and find indices associated with the mt_array. 
                idxA = find( ~cellfun( @isempty, strfind( vecLabels, strAster ) ) );
                vecA = vec( idxA );
                vecLabelsA = erase( vecLabels( idxA ), strAster );
                obj.featureList{ jAster} = absorbVec( obj.featureList{ jAster}, vecA, vecLabelsA );

            end
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            % Order of features : 
            % AsterMTOC (complete)
            % Aster MTs (poly coef end)

            vecList = {};
            vecLabelsList = {};
            objList = {};

            % iAster vectors
            nFeat = 0;
            for jAster = 1 : obj.numFeatures

                nFeat=nFeat+1;
                % get mtoc vector
                [ vecList{nFeat} , vecLabelsList{nFeat}] = getVec( obj.featureList{jAster}.featureList{1}, obj.props2Fit.fit{obj.dim}.aster.spot);
                objList{nFeat} = obj.featureList{jAster}.featureList{1};

                % get imt vectors
                nummt = length( obj.featureList{jAster}.featureList) - 1;
                for jmt = 1 : nummt
                    [ vecList{ nFeat+jmt} , vecLabelsList{ nFeat+jmt}] = getVec( obj.featureList{jAster}.featureList{1+jmt}, obj.props2Fit.fit{obj.dim}.aster.curve);
                    objList{ nFeat+jmt} = obj.featureList{jAster}.featureList{1+jmt};
                end
                nFeat = length( vecList);
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
            
            % In an IMT bank with N asters, we can have, at maximum, Nx2
            % curve features. For each aster, we will find a possible microtubule to add.
            % Then we will pick the best of the possibilities and add them to our feature list.

            successAdd = 0;
            polyOrder = [2 2 1];
            try
                polyOrder = obj.polyOrder;
            end

            if obj.dim == 2
                sigma=[1.2 1.2];
            elseif obj.dim == 3
                sigma=[1.2 1.2 1.0];
            end

            % Ask existing asters to find missing features
            feature = cell( 1, obj.numFeatures);
            residual_density = [];
            for jFeature = 1 : obj.numFeatures
                
                nCurves = length( obj.featureList{jFeature}.featureList)-1;
                
                switch nCurves
                    case 2
                        continue
                    case 0
                        error('An aster should never have 0 curves')
                    case 1
                        
                        % Make image without all other asters
                        imElse = 0*Image2Find;
                        for jFeature2 = 1 : obj.numFeatures
                            if jFeature ~= jFeature2
                                imElse = imElse + obj.featureList{ jFeature2}.simulateFeature( size(Image2Find));
                            end
                        end
                        imElse = imElse < 0.1*max(imElse(:));
                        
                        % Find start point of missing curve: its the
                        % position of the mtoc
                        stPoint = obj.featureList{ jFeature}.featureList{1}.position;
                        
                        % Find curve coords of given orientation. enforce
                        % having atleast 2 coordinates
                        coords = Curve.findCurve( Image2Find.*imElse, stPoint(1:2), obj.featureList{ jFeature}.featureList{2}.GetOrientation()+pi );
                        if size( coords, 2) < 2
                            continue
                        end
                        if obj.dim == 3
                            coords(3,:) = stPoint(3);
                        end
                        
                        % Find coefficients
                        coeffs = Curve.estimatePolyCoefficients( coords, polyOrder);
                        if obj.dim == 3 && all( coords(3,:) == size(Image2Find,3) )
                            coeffs{3} = -coeffs{3};
                        end

                        % Find amplitudes
                        bkg = median( obj.image(:) );
                        curveAmp = median( Cell.findAmplitudeAlongCurve( max( obj.image, [], 3), coeffs ) ) - bkg;
                        if curveAmp < 0
                            curveAmp = bkg;
                            warning( 'curveAmp was less than bkg')
                        end

                        % Create curve features
                        feature{jFeature} = Curve( stPoint, coeffs, curveAmp, sigma, obj.dim, obj.props2Fit.fit{obj.dim}.curve, ...
                            obj.props2Fit.graphics.curve);
                        feature{jFeature}.findVoxelsInsideMask( logical(obj.image) );
                        
                        % Find length and unit residual of prospective feature
                        len = feature{jFeature}.GetLength();
                        imFeat2D = max( feature{jFeature}.simulateFeature( size(obj.image) ), [], 3);
                        imFeat2DMask = imFeat2D > 0.1*max(imFeat2D(:));
                        Image2Find2D = max(Image2Find,[],3);
                        residual_density(jFeature) = sum( Image2Find2D(:) .* imFeat2DMask(:) ) / len;
                        
                end
            end
            % find the best feature to add
            [~, idx] = max( residual_density);
            
            % Add the best feature
            if ~isempty(idx)
                idxAdd = 1+ length( obj.featureList{idx}.featureList);
                obj.featureList{idx}.addFeatureToList( feature{idx} );
                % Add feature ID and update the featureMap
                feature{idx}.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
                obj.featureMap( feature{idx}.ID) = [ idx idxAdd];
                successAdd = 1;
            end
                            
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
            
            % Ask asters to give their worst microtubule features
            feat = cell( 1, obj.numFeatures);
            succ = zeros( 1, obj.numFeatures);
            for jFeature = 1 : obj.numFeatures
                [feat{jFeature}, succ(jFeature) ] = obj.featureList{jFeature}.findLeastPromFeature();
            end
            
            % If no features to remove, exit the function
            if all( succ == 0)
                successRemove = 0;
                return
            end

            % Find least prominent feature
            p = 1e6; idxAster = 0; idx = 0;
            for jFeature = 1 : obj.numFeatures
                if succ(jFeature) && feat{jFeature}.p < p
                    p = feat{jFeature}.p;
                    idxAster = jFeature;
                    idx = feat{jFeature}.idx;
                end
            end          

            % Remove the worst microtubule
            obj.featureList{ idxAster}.removeFeatureFromList( idx);
            
            % If, as a result, an aster contains only a spot, remove it
            % with the microtubule
            if obj.featureList{ idxAster}.numFeatures == 1
                obj.removeFeatureFromList( idxAster);
            end
                
        end
        % }}}
        
        % addOrganizers {{{
        function [obj, successAdd] = addOrganizers( obj, Image2Find)
            % Searches for missing subfeatures and adds them to the featureList
            
            % In an IMT bank with N asters, we can have, at maximum, Nx2
            % curve features. For each aster, we will find a possible microtubule to add.
            % Then we will pick the best of the possibilities and add them to our feature list.

            polyOrder = [2 2 1];
            try
                polyOrder = obj.polyOrder;
            end

            minLength = 10;
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
            
            % If useful curve found, find the Seed/MTOC, find curves and form an aster.
            % Find Seed/MTOC
            [cX,cY] = Methods.InterpolateCoords( coord(1,:), coord(2,:), 20); 
            indXY = sub2ind( size(Image2Find2D), round(cY), round(cX) );
            [IntSeed, IdxSeed] = max( Image2Find2D( indXY) );
            [ Seed(2), Seed(1)] = ind2sub( size(Image2Find2D), indXY(IdxSeed) );
            
            % Return if curve seed is not bright
            if max( obj.image( Seed(2), Seed(1), :) ) < stats.ThreshLow
                successAdd = 0;
                return
            end
            
            % add 3D coordinate if dimension is 3D
            if obj.dim == 3
                Seed(3) = idxZ( Seed(2), Seed(1) );
%                 coord(3,:) = Seed(3);
            end
            mtoc = Spot( Seed, IntSeed, sigma, obj.dim, obj.props2Fit.fit{obj.dim}.spot, ...
                obj.props2Fit.graphics.spot);
            
            % Break curve into 2 mts
            mts{1} = [ cX( IdxSeed:-1:1); cY( IdxSeed:-1:1) ];
            mts{2} = [ cX( IdxSeed:end); cY( IdxSeed:end) ];
            if obj.dim == 3
                mts{1}(3,:) = Seed(3);
                mts{2}(3,:) = Seed(3);
            end
            
            % Create interphase MTs
            iMT = cell( 1, length( mts ) );
            for jmt = 1 : length( mts )

                % Get Coefficients of iMT
                coeffs = Curve.estimatePolyCoefficients( mts{jmt}, polyOrder);
                
                if obj.dim == 3 && all( mts{jmt}(3,:) == size(Image2Find,3) )
                    coeffs{3} = -coeffs{3};
                end
                curveAmp = median( Cell.findAmplitudeAlongCurve( Image2Find2D, coeffs ) );
                if curveAmp < 0
                    curveAmp = 0.05;
                    warning( 'curveAmp was less than 0')
                end

                iMT{jmt} = Curve( mtoc.position, coeffs, curveAmp, sigma, obj.dim, obj.props2Fit.fit{obj.dim}.curve, ...
                    obj.props2Fit.graphics.curve);
                iMT{jmt}.findVoxelsInsideMask( logical(Image2Find) );

            end

            % Store iMTOC + iMTs in Aster Objects
            aster = AsterMT( obj.dim, mtoc, iMT{:} );
            
            % Add the aster
            obj.addFeatureToList( aster );
            % Add feature ID and update the featureMap
            aster.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            aster.syncFeatures();
            successAdd = 1;
                            
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
        
        % getSubFeatureNumber {{{
%         function numFeat = getSubFeatureNumber( obj)
%             numFeat = length( obj.featureList);
%         end
        % }}}
        
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
                        feats( 
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
            nAsters = length( coords);
            
            % Create the MTOCs
            MTOC = cell(1, nAsters); 
            [mip, ~] = max( imageIn, [],3);
            for jmtoc = 1 : nAsters
                coord = coords{jmtoc}{1}(:,1);
                amp = mip( coord(2), coord(1) ) - bkg;
                MTOC{ jmtoc} = Spot( coord', amp, sigma, dim, ...
                    props.fit{dim}.spot, props.graphics.spot);
            end

            % Create interphase MT and interphase Aster objects
            iAsters = cell(1, nAsters );
            for jMTOC = 1 : nAsters

                nMT = length( coords{jMTOC} );
                % Create interphase MTs
                iMT = cell( 1, nMT); 
                for jmt = 1 : nMT
                    
                    % Get Coefficients of microtubules
                    coeffs = Curve.estimatePolyCoefficients( coords{jMTOC}{jmt}, params.PolyOrder);
                    
                    % Ensure coefficients are within the image region
                    if all( coords{jMTOC}{jmt}(3,:) == size(imageIn,3) )
                        coeffs{3}(1:end-1) = -abs(coeffs{3}(1:end-1));
                    end
                    
                    % Ensure coefficients are within the image region
                    if all( coords{jMTOC}{jmt}(3,:) == 1 )
                        coeffs{3}(1:end-1) = +abs(coeffs{3}(1:end-1));
                    end
                    
                    % Get amplitude of microtubule
                    curveAmp = median( Cell.findAmplitudeAlongCurve( mip, coeffs ) ) - bkg;
                    if curveAmp < 0
                        curveAmp = bkg;
                        warning( 'IntBank.findIntBank: curveAmp was less than bkg. Forced equality.')
                    end
                    
                    % Create 
                    iMT{jmt} = Curve( MTOC{jMTOC}.position, coeffs, curveAmp, sigma, dim, props.fit{dim}.curve, props.graphics.curve);
                    iMT{jmt}.findVoxelsInsideMask( logical(imageIn) );

                end

                % Store iMTOC + iMTs in Aster Objects
                iAsters{ jMTOC} = AsterMT( dim, MTOC{ jMTOC}, iMT{:} );

            end

            % Create the interphase aster bank.
            intBankObj = IMTBank( dim, imageIn, iAsters, props);
            intBankObj.polyOrder = params.PolyOrder;
            intBankObj.findEnvironmentalConditions();
            intBankObj.forceInsideMask();

        end
        % }}}

    end
end
