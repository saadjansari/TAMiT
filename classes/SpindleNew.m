classdef SpindleNew < OrganizerMaster
% A Spindle is a higher level feature that sits inside a mitotic cell. It is composed of 2 MT_arrays and a Microtubule connecting the 2 arrays
    properties
    end

    methods (Access = public)
        
        % Spindle {{{
        function obj = SpindleNew( dim, image, featureList, props2Fit)

            if ~strcmp( featureList{1}.type, 'Line')
                error('Spindle: featureList{1} must be of type ''Line'' ') 
            end

            obj = obj@OrganizerMaster( dim, image, featureList, props2Fit, 'SpindleNew');
            obj.numFeatures = length(featureList);
            obj.mask = ones( size(obj.image));

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels, ub, lb] = getVec( obj, props)

            vec = []; lb = []; ub = []; vecLabels={};
            % Get general environmental properties
            [vecE, vecLabelsE] = obj.getVecEnvironment();
            ubE = vecE; lbE = vecE;

            % get vectors from features
            props.spindle = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            props.spb = {'position', 'amplitude', 'sigma'};
            props.mt = {'thetaInit', 'normalVec', 'L', 'amplitude', 'sigma'};
            props.mt2 = {'thetaInit', 'curvature', 'L', 'amplitude', 'sigma'};
            warning('using hard coded curvature property') 
            
            % get vector from Spindle microtubule
            [vec_spindle, vecLabels_spindle, ub_s, lb_s] = getVec( obj.featureList{ 1}, props.spindle );
            vec = [vec, vec_spindle]; ub = [ub, ub_s]; lb = [lb, lb_s];
            vecLabels_spindle = strcat( 'S_', vecLabels_spindle);
            vecLabels = { vecLabels{:}, vecLabels_spindle{:} };

            % Loop over mt arrays and get their vectors. Remove the SPBs from the fit vectors. 
            for jmt = 2 : obj.numFeatures
                [vec_mtarray, vecLabels_mtarray, ub_ma, lb_ma] = getVec( obj.featureList{ jmt}, props.spb, props.mt2);
                % find indexes of strings that start with the 'SPB_position'
                idxRm = find( strcmp( vecLabels_mtarray, 'SPB_position') );
                vec_mtarray( idxRm) = []; ub_ma( idxRm) = []; lb_ma( idxRm) = [];
                vecLabels_mtarray( idxRm) = [];
                vec = [vec , vec_mtarray]; ub = [ub, ub_ma]; lb = [lb, lb_ma];
                vecLabels_mtarray = strcat( 'SP_', num2str( jmt-1), '_', vecLabels_mtarray);
                vecLabels = { vecLabels{:}, vecLabels_mtarray{:} };
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels, errBoolean)
            
            if nargin < 4
                errBoolean = 0;
            end
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels, errBoolean);

            % Spindle Vector
            % Take the vector and find the indexes associated with the spindle parameters
            idxSpindle = find( ~cellfun( @isempty, strfind( vecLabels, 'S_') ) );
            idxSpindlePosition( 1, :) = find( ~cellfun( @isempty, strfind( vecLabels, 'S_startPosition') ) );
            idxSpindlePosition( 2, :) = find( ~cellfun( @isempty, strfind( vecLabels, 'S_endPosition') ) );

            % Get the vector for the spindle by removing the spindle substring. Absorb the vector
            vecS = vec( idxSpindle);
            vecLabelsS = erase( vecLabels( idxSpindle), 'S_');
            obj.featureList{ 1} = absorbVec( obj.featureList{ 1}, vecS, vecLabelsS, errBoolean );

            % Aster Vector
            % Create the vector for each subfeature MT_Array and send it to the objects for absorption
            for jAster = 1 : obj.numFeatures-1

                strAster = [ 'SP_', num2str(jAster), '_' ];
                % Take the vector and find indices associated with the mt_array. 
                idxSP = find( ~cellfun( @isempty, strfind( vecLabels, strAster ) ) );
                vecSP = vec( idxSP );
                vecLabelsSP = erase( vecLabels( idxSP ), strAster );
                % Append the parameters for the Spindle Pole Body for this array
                vecSP = [ vec( idxSpindlePosition( jAster, :) ) , vecSP ];
                vecLabelsSPB = repmat( {'SPB_position'}, 1, length( idxSpindlePosition( jAster, :) ) );
                vecLabelsSP = { vecLabelsSPB{:}, vecLabelsSP{:} };
                obj.featureList{ 1+jAster} = absorbVec( obj.featureList{ 1+jAster}, vecSP, vecLabelsSP, errBoolean );

            end

        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList, ubList, lbList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            % Order of features : 
            %   Spindle Microtubule (complete optimization)
            %   Astral Microtubules (end optimization)

            vecList = {};
            vecLabelsList = {};
            objList = {};
            ubList = {}; lbList = {};

            % Define properties to get
            props.spindle = {'startPosition', 'endPosition', 'sigma'};
            props.spb = {'sigma'};
            props.mt = {'thetaInit', 'normalVec', 'L', 'amplitude', 'sigma'};

            % Spindle vector
            [ vecList{1} , vecLabelsList{1}, ubList{1}, lbList{1}] = getVec( obj.featureList{1}, props.spindle );
            objList{1} = obj.featureList{1};
            numFeatures = 1;

%             % SPB 1 vector
%             [ vecList{2} , vecLabelsList{2}, ubList{2}, lbList{2}] = getVec( obj.featureList{2}.featureList{1}, props.spb);
%             objList{2} = obj.featureList{2}.featureList{1};
%             numFeatures = numFeatures+1;
% 
%             % SPB 2 vector
%             [ vecList{3} , vecLabelsList{3}, ubList{3}, lbList{3}] = getVec( obj.featureList{3}.featureList{1}, props.spb);
%             objList{3} = obj.featureList{3}.featureList{1};
%             numFeatures = numFeatures+1;

            % Microtubule vector : Aster 1
            for jmt = 1 : obj.featureList{2}.numFeatures-1
                [ vecList{ numFeatures+jmt} , vecLabelsList{ numFeatures+jmt}, ubList{numFeatures+jmt}, lbList{numFeatures+jmt}] = getVec( obj.featureList{2}.featureList{1+jmt}, props.mt );
                objList{ numFeatures+jmt} = obj.featureList{2}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{2}.numFeatures-1;

            % Microtubule vector : Aster 2
            for jmt = 1 : obj.featureList{3}.numFeatures-1
                [ vecList{ numFeatures+jmt} , vecLabelsList{ numFeatures+jmt}, ubList{numFeatures+jmt}, lbList{numFeatures+jmt}] = getVec( obj.featureList{3}.featureList{1+jmt}, props.mt );
                objList{ numFeatures+jmt} = obj.featureList{3}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{3}.numFeatures-1;

%             % Prepend environmental parameters to each vec and vecLabel
%             [vecE, vecLabelsE] = obj.getVecEnvironment();
%             
%             for jFeat = 1 : numFeatures
%                 vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
%                 vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
%             end

        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures( obj)

            % Get the spindle microtubule start and end positions and update the SPB associated with the spindle
            obj.featureList{2}.featureList{1}.position = obj.featureList{1}.startPosition;
            obj.featureList{3}.featureList{1}.position = obj.featureList{1}.endPosition;

        end
        % }}}

        % getAngleXY {{{
        function angle = getAngleXY( obj)
            % Return the angle the spindle makes in XY : 2 values correspond to the two origins at the two poles

            angle(1) = mod( atan2( obj.featureList{1}.endPosition(2) - obj.featureList{1}.startPosition(2) , obj.featureList{1}.endPosition(1) - obj.featureList{1}.startPosition(1) ), 2*pi );

            angle(2) = mod( atan2( obj.featureList{1}.startPosition(2) - obj.featureList{1}.endPosition(2) , obj.featureList{1}.startPosition(1) - obj.featureList{1}.endPosition(1) ), 2*pi );

        end
        % }}}

        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, Image2Find)
            % Searches for missing features and adds them to the featureList
            
            % For a spindle, we can only add more astral microtubules
            % For each pole we will find a possible microtubule to add, then we will pick the best of the two possibilities and add them to our feature list.

            % How to find a possible microtubule:
            %   1. Angular Sweep (not clear how to ensure a line is found if there are no peaks)
            %   2. Pick the location of the highest residual and draw a line from both pole connecting to the max residual point. Find the mean intensity of each line. Subtract the mean from the actual voxel values along the line to find a total residual. Pick the pole with the minimum residual.

            % Get Spindle Angle
            spindleAngle = getAngleXY( obj);
            spindleExclusionRange = deg2rad( 60);

            % ask subfeatures to find missing features
            [feature{1}, suc1, res1] = obj.featureList{2}.findMissingFeature( obj.image, Image2Find, spindleAngle(1), spindleExclusionRange );
            [feature{2}, suc2, res2] = obj.featureList{3}.findMissingFeature( obj.image, Image2Find, spindleAngle(2), spindleExclusionRange );

            % Make higher level decision on the best missing feature
            if suc1 && ~suc2
                idxKeep = 1;
            elseif suc2 && ~suc1
                idxKeep = 2;
            elseif ~suc1 && ~suc2
                successAdd = 0; return
            else
                [~,idxKeep] = min( [res1, res2]);
            end

            % Add feature to the correct subfeature 
            idxAdd = obj.featureList{1+idxKeep}.addFeatureToList( feature{idxKeep} );
            imbkg = 0*obj.image; imbkg(:) = obj.background;
            feature{idxKeep}.preOptimize(obj.image, imbkg);

            % Add feature ID and update the featureMap
            feature{idxKeep}.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            obj.featureMap( feature{idxKeep}.ID) = [ 1+idxKeep idxAdd];
            successAdd = 1;
        end
        % }}}

        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % For a spindle, we can only remove  astral microtubules
            % For each pole we will find a possible microtubule to remove, then we will pick the best of the two possibilities and remove them from our feature list.
            successRemove = 1;
            % Ask asters to give their worst microtubule features
            [feat{1}, succ1] = obj.featureList{2}.findHighResFeature(Image2Find);
            [feat{2}, succ2] = obj.featureList{3}.findHighResFeature(Image2Find);

            % If no features to remove, exit the function
            if succ1 == 0 && succ2 == 0
                successRemove = 0;
                return
            % Make higher level decision on the worst-est feature, find least prominent feature
            elseif succ1 == 0
                idxAster=2;
            elseif succ2==0
                idxAster=1;
            else
                [~, idxAster ] = max( [ feat{1}.residual, feat{2}.residual] );
            end

            % Remove the worst microtubule
            obj.featureList{ 1+idxAster }.removeFeatureFromList( feat{idxAster}.idx);
                
        end
        % }}}

        % finalizeAddedFeatures {{{
        function obj = finalizeAddedFeatures(obj)
            
            % Display properties for lines in asters
            display = {'Color', [1 0 1], 'LineWidth', 3};
            
            if length( obj.featureList ) > 1
                for jL = 2 : length(obj.featureList{2}.featureList )
                    obj.featureList{2}.featureList{jL}.display = display;
                end
                for jL = 2 : length(obj.featureList{3}.featureList )
                    obj.featureList{3}.featureList{jL}.display = display;
                end
            end
        end
        % }}}
        
        % findEnvironmentalConditions {{{
        function obj = findEnvironmentalConditions( obj)
            % Finds the cytoplasmic and the nucleoplasmic backgrounds

            env = obj.props2Fit.fit{obj.dim}.Environment;

            % Check which environmental conditions exist
            isBkg = any( cellfun( @(x) strcmp(x,'background'), env ));
            isBkgNuc = any( cellfun( @(x) strcmp(x,'backgroundNuclear'), env ) );

            % Background
            if isBkg
                obj.mask = logical(obj.image);
                imVals = obj.image( find(obj.mask));
                obj.background = median( imVals);
            end

            % Nuclear Background
            if isBkgNuc
                
                % Nuclear mask
                mask = BY_find_nuclear_mask( obj.image);
                bkg_nuc = median( obj.image( find(mask(:) ) ) );
                
                % new background first
                if isBkg
                    obj.background = median( obj.image( obj.mask & ~mask ) );
                end
                obj.backgroundNuclear = bkg_nuc - obj.background;
                obj.maskNuclear = mask;
            end

        end
        % }}}
        
    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'SpindleNew')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = Line.loadFromStruct( S.featureList{ 1} ); 
            for jFeat = 2 : length( S.featureList)
                featureList{ jFeat} = Aster.loadFromStruct( S.featureList{ jFeat} ); 
            end

            obj = SpindleNew( S.dim, S.image, featureList, S.props2Fit);
%             obj = obj@Organizer( S.dim, featureList, S.type);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap();

        end
        % }}}

        % findSpindle {{{
        function spindleObj = findSpindle( imageIn, params, props)
            % Finds a spindle in the given image

            % Get dimensionality and background
            dim = length( size(imageIn) );
            bkg = median( imageIn( imageIn(:) > 0) );

            % Define sigma for this bank
            if dim == 2
                sigma=[1.2 1.2];
            elseif dim == 3
                sigma=[1.2 1.2 1.0];
            end
            spindleObj = 1;
        end
        % }}}
    end

end
