classdef Aster < Organizer
% This is a curved Aster. This sits inside a spindle. It contains 1 SPB (Spindle Pole Body) and a number of curved microtubules
    properties
    end

    methods
       
        % Aster {{{
        function obj = Aster( dim, spb, varargin)

            if nargin < 2 
                error( 'MT_Array: number of inputs is inconsistent')
            end

            obj = obj@Organizer( dim, { spb, varargin{:} }, 'Aster');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels,ub,lb] = getVec( obj, propsSPB, propsMT)

            if nargin < 3
                propsMT = {'thetaInit', 'curvature','L','amplitude', 'sigma'};
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

                [vec_mt, vecLabels_mt, ub_mt, lb_mt] = getVec( obj.featureList{ jmt}, propsMT );
                vec = [vec , vec_mt]; lb = [lb, lb_mt]; ub = [ub, ub_mt];
                vecLabels_mt = strcat( 'MT', num2str( jmt-1), '_', vecLabels_mt);
                vecLabels = { vecLabels{:}, vecLabels_mt{:} };
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

            % Get the vector for each microtubule, and ask the microtubule objects to absorb the vector
            for jmt = 2 : obj.numFeatures
                strMT = ['MT', num2str(jmt-1), '_'];
                idxMT = find( ~cellfun( @isempty, strfind( vecLabels, strMT ) ) );
                vecMT = vec( idxMT );
                vecLabelsMT = erase( vecLabels( idxMT), strMT );
                obj.featureList{ jmt} = absorbVec( obj.featureList{ jmt}, vecMT, vecLabelsMT, errBoolean);
            end
            
            p0 = obj.featureList{ 1}.position;
            obj.featureList{ 1} = absorbVec( obj.featureList{ 1}, vecSPB, vecLabelsSPB, errBoolean );
            p1 = obj.featureList{ 1}.position;
            
            for jmt = 2 : obj.numFeatures
                obj.featureList{jmt}.origin = obj.featureList{jmt}.origin + (p1-p0);
            end

        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures( obj )

            % Get the SPB  position and update the MTs  associated with the SPB
            for jmt = 1 : obj.numFeatures-1
%                 obj.featureList{1+jmt}.origin = obj.featureList{1}.position;
                obj.featureList{1+jmt}.origin = obj.featureList{1}.position + 5*[cos(spindleAngle(1)+pi), sin(spindleAngle(1)+pi)];
            end

        end
        % }}}

        % findMissingFeature {{{
        function [mt, success, resid] = findMissingFeature( obj, imOrg, Image2Find, xAngle, xRange, featType)

%             props = {'origin', 'thetaInit', 'normalVec','L', 'amplitude', 'sigma'};
            props = {'origin', 'thetaInit', 'normalVec','L', 'amplitude', 'sigma'};
            display = {'Color', [1 0 0], 'LineWidth', 3};
            if obj.dim==3, sigma=[2.5 2.5 1]; elseif obj.dim==2, sigma=[2.5 2.5]; end
            Image2Find( Image2Find < 0) = 0; 
            im3g = imgaussfilt3(Image2Find,1);
            im2g = max(im3g,[],3);
                                    
            minLength = 10;
            [Image2Find2D, idxZ] = max(Image2Find, [],3);
            
            % Do a sweep radial sweep and look for peaks corresponding to
            % rods.            
            startPoint = obj.featureList{1}.position(1:2) + 3*[ cos(xAngle+pi), sin(xAngle+pi)];
            [phiIntensity, phiValues] = Cell.radIntegrate2D( im2g, startPoint, 5, 25);
            
            % find peaks in Phi Intensity
            nbkg = median(im2g(:))*2;
            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, ...
                'MinPeakHeight', 2*nbkg, 'MinPeakProminence', nbkg);
            peakPhi = mod( peakPhi, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )
            
            nms3 = 0*imOrg;
            for jZ = 1 : size(imOrg,3)
                [~, ~, nms3(:,:,jZ), ~] = steerableDetector( imgaussfilt(imOrg(:,:,jZ),1), 4, 3);
            end
            
            % find peaks in Phi Intensity
            if length(peakPhi) == 0
                success = 0; mt = []; resid=[];
                return
            end

            % remove angle belonging to spindle.
            spindleAngle = mod( xAngle, 2*pi); 
            idxRm = find( abs(peakPhi-spindleAngle) < xRange | abs(peakPhi-spindleAngle+2*pi) < xRange | abs(peakPhi-spindleAngle-2*pi) < xRange);
            peakPhi( idxRm)=[];

            if length(peakPhi) == 0
                success = 0; mt = []; resid=[];
                return
            end
            
            % Remove angle belonging to pre-existing curves
            idxRm = [];
            for j1 = 1 : length(peakPhi)
                for j2 = 2 : obj.numFeatures
                    phi1 = peakPhi(j1);
                    % get phi of pre-existing line
                    cc2 = obj.featureList{j2}.GetCoords();
                    phi2 = atan2( cc2(2,2)-cc2(2,1) , cc2(1,2)-cc2(1,1) );
                    if abs( phi1 - phi2) < 0.3 || abs( phi1 - phi2 -2*pi) < 0.3 || abs( phi1 - phi2 + 2*pi) < 0.3
                        idxRm = [idxRm, j1];
                    end
                end
            end
            peakPhi( idxRm)=[];

            if length(peakPhi) == 0
                success = 0; mt = []; resid=[];
                return
            end
            
            % Default Values
            defaultVisibility = 10;
            defaultFieldOfView = 70;
            defaultStepSize = 4;
            minLength = 10;
            
            nms3 = 0*Image2Find;
            for jZ = 1 : size(Image2Find,3)
                [~, ~, nms3(:,:,jZ), ~] = steerableDetector(Image2Find(:,:,jZ), 4, 2);
            end
            
            % create a mask to apply to steerable image
            mask = im2g;
            st = Methods.GetImageStats(im2g,0);
            mask( mask < 2*st.Median) = 0; mask(mask ~=0) = 1;
            
            % Find curved MT 
            missingFeatures = {};
            for pp = 1 : length(peakPhi)
                               
                % Set the two opposite angles for search
                cc = Methods.estimateCurveCoords( startPoint', peakPhi(pp), sum(nms3,3).*max(mask,[],3), defaultStepSize, defaultVisibility, defaultFieldOfView, 1);
                
                % figure out z-coordinates
                [a2D, i2d] = max( imgaussfilt(Image2Find, 1), [], 3);
                c3 = zeros( 1, size(cc,2));
                for jj = 1 : length(c3)
                    rng = 1; alist = []; idlist = [];
                    for j1 = -rng : rng
                        for j2 = -rng : rng
                            try
                                alist = [ alist, a2D( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                                idlist = [ idlist, i2d( round(cc(2,jj))+j1, round(cc(1,jj))+j2)];
                            end
                        end
                    end
                    [~, idd] = max( alist); c3( jj) = idlist(idd);
                end
                cc = [cc; c3];
                
                % Length of end-end curve
                if norm( cc(:,end)-cc(:,1)) > minLength
                    missingFeatures = { missingFeatures{:}, cc};
                end
                
            end
            
            % Remove any duplicates
            phiz = zeros(1, length(missingFeatures));
            for j1 = 1 : length( missingFeatures)
                phiz(j1) = atan2( missingFeatures{j1}(2,2)-missingFeatures{j1}(2,1) , missingFeatures{j1}(1,2)-missingFeatures{j1}(1,1) );
            end
            destroy = [];
            for p1 = 1:length(phiz)
                for p2 = 1:length(phiz)
                    if p1 ~=p2
                        if abs( phiz(p1) - phiz(p2)) < 0.1 || abs( phiz(p1) - phiz(p2) -2*pi) < 0.1 || abs( phiz(p1) - phiz(p2) + 2*pi) < 0.1
                            % figure out which one to destroy
                            if norm( missingFeatures{p1}(:,end)-missingFeatures{p1}(:,1)) > norm( missingFeatures{p2}(:,end)-missingFeatures{p2}(:,1))
                                destroy = [destroy, p2];
                            else
                                destroy = [destroy, p1];
                            end
                        end
                    end
                end
            end
            destroy = unique(destroy);
            missingFeatures( destroy) = [];
            

            % Create curved MT objects
            curvedMTs = cell(1, length( missingFeatures ));
            idxRm = [];
            for jb = 1 : length( curvedMTs )

                coords = missingFeatures{jb};
                        
                % get intial theta, curvature coefficients from
                % coords
                curvature_model = obj.parameters.curvature_model;
                [origin,thetaInit,curvature, L, ccd] = get_tan_curvature_coeffs( coords, curvature_model);

                % Ensure origin is within the Z stacks
                if origin(3) >= size( Image2Find,3)
                    origin(3) = size(Image2Find,3)-0.2;
                elseif origin(3) <= 1
                    origin(3) = 1.2;
                end

                % Get amplitude
                A1 = smooth( findAmplitudeAlongCurveCoords( max(Image2Find,[],3), round(ccd(1:2,:)) ));
                A1( A1 < 0) = 0; amp = median(A1);

                % Ensure initial theta is pointed to inside the image region
                if obj.dim == 3
                    if origin(3) == 1
                        thetaInit(1,2) = pi/2 - 0.03;
                    elseif origin(3) == size(Image2Find, 3)
                        thetaInit(1,2) = pi/2 + 0.03;
                    end
                end

                curvedMTs{jb} = CurvedMT( curvature_model, origin', thetaInit, curvature, L, amp, sigma, obj.dim, props, display);
                    
            end
            curvedMTs(idxRm) = [];
                
            % remove curves shorter than minLength
            idxRm = [];
            for jb = 1 : length(curvedMTs)
                if curvedMTs{jb}.GetLength() < minLength
                    idxRm = [idxRm, jb];
                end
            end
            curvedMTs( idxRm) = [];

            if length( curvedMTs) == 0
                success = 0; mt = []; resid=[];
                return
            end

            % find best curved microtubule
            res = [];
            for jb = 1 : length(curvedMTs)
                imSim = curvedMTs{jb}.simulateFeature(size(Image2Find));
                res = [res, sum( (Image2Find(:) -imSim(:) ).^2 )];
            end
            [resid, idxSelect] = min( res);
            mt = curvedMTs{ idxSelect};
            mt.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            success = 1;

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
                res(jmt) = sum( (Cell.findAmplitudeAlongCurveCoords( ImageRef, round( obj.featureList{jmt}.GetCoords()))).^2 );
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
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Aster')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = Spot.loadFromStruct( S.featureList{ 1} ); 
            for jFeat = 2: length( S.featureList)
                try
                    featureList{ jFeat} = Line.loadFromStruct( S.featureList{jFeat} );
                catch
                    featureList{ jFeat} = CurvedMT.loadFromStruct( S.featureList{jFeat} );
                end
            end

            obj = Aster( S.dim, featureList{:});

        end
        % }}}
        
    end

end
