classdef Spindle < Organizer
% A Spindle is a higher level feature that sits inside a mitotic cell. It is composed of 2 MT_arrays and a Microtubule connecting the 2 arrays
    properties
        numAsters % numAsters = 2
        props2Fit
        background
        backgroundNuclear 
        maskNuclear
        featureMap
        image
    end

    methods (Access = public)
        
        % Spindle {{{
        function obj = Spindle( dim, image, featureList, props2Fit)

            % ensure that featureList has 3 elements. The first is of type Line, the second and third are of type AsterMT
            if length( featureList) ~=3
                error('Spindle: featureList input argument must contain 3 objects in total')
            end
            if ~strcmp( featureList{1}.type, 'Line')
                error('Spindle: featureList{1} must be of type ''Line'' ') 
            end
            if ~strcmp( featureList{2}.type, 'AsterMT')
                error('Spindle: featureList{2} must be of type ''AsterMT'' ') 
            end
            if ~strcmp( featureList{3}.type, 'AsterMT')
                error('Spindle: featureList{3} must be of type ''AsterMT'' ') 
            end

            obj = obj@Organizer( dim, featureList, 'Spindle');
            obj.image = image;
            obj.props2Fit = props2Fit;
            
            % initialize featureMap and assign IDs
            obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            % assign voxels to its basic elements
            obj.featureList{1}.findVoxelsInsideMask( logical(image) );
            for jAster = 1 : 2
                for jFeat = 1 : length( obj.featureList{ 1+jAster}.featureList)

                    obj.featureList{1+jAster}.featureList{jFeat}.findVoxelsInsideMask( logical( image) );

                end
            end

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props)

            % Get general environmental properties
            [vec, vecLabels] = getVecEnvironment( obj, obj.props2Fit);

            % get vectors from features
            props.spindle = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            props.spb = {'position', 'amplitude', 'sigma'};
            props.mt = {'endPosition', 'amplitude', 'sigma'};

            % get vector from Spindle microtubule
            [vec_spindle, vecLabels_spindle] = getVec( obj.featureList{ 1}, props.spindle );
            vec = [vec, vec_spindle];
            vecLabels_spindle = strcat( 'S_', vecLabels_spindle);
            vecLabels = { vecLabels{:}, vecLabels_spindle{:} };

            % Loop over mt arrays and get their vectors. Remove the SPBs from the fit vectors. 
            for jmt = 2 : obj.numFeatures
                [vec_mtarray, vecLabels_mtarray] = getVec( obj.featureList{ jmt}, props.spb, props.mt);
                % find indexes of strings that start with the 'SPB_position'
                idxRm = find( strcmp( vecLabels_mtarray, 'SPB_position') );
                vec_mtarray( idxRm) = [];
                vecLabels_mtarray( idxRm) = [];
                vec = [vec , vec_mtarray];
                vecLabels_mtarray = strcat( 'SP_', num2str( jmt-1), '_', vecLabels_mtarray);
                vecLabels = { vecLabels{:}, vecLabels_mtarray{:} };
            end

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels);

            % Spindle Vector
            % Take the vector and find the indexes associated with the spindle parameters
            idxSpindle = find( ~cellfun( @isempty, strfind( vecLabels, 'S_') ) );
            idxSpindlePosition( 1, :) = find( ~cellfun( @isempty, strfind( vecLabels, 'S_startPosition') ) );
            idxSpindlePosition( 2, :) = find( ~cellfun( @isempty, strfind( vecLabels, 'S_endPosition') ) );
            % Get the vector for the spindle by removing the spindle substring. Absorb the vector
            vecS = vec( idxSpindle);
            vecLabelsS = erase( vecLabels( idxSpindle), 'S_');
            obj.featureList{ 1} = absorbVec( obj.featureList{ 1}, vecS, vecLabelsS );
             
            % Aster Vector
            % Create the vector for each subfeature MT_Array and send it to the objects for absorption
            for jAster = 1 : 2

                strAster = [ 'SP_', num2str(jAster), '_' ];
                % Take the vector and find indices associated with the mt_array. 
                idxSP = find( ~cellfun( @isempty, strfind( vecLabels, strAster ) ) );
                vecSP = vec( idxSP );
                vecLabelsSP = erase( vecLabels( idxSP ), strAster );
                % Append the parameters for the Spindle Pole Body for this array
                vecSP = [ vec( idxSpindlePosition( jAster, :) ) , vecSP ];
                vecLabelsSPB = repmat( {'SPB_position'}, 1, length( idxSpindlePosition( jAster, :) ) );
                vecLabelsSP = { vecLabelsSPB{:}, vecLabelsSP{:} };
                obj.featureList{ 1+jAster} = absorbVec( obj.featureList{ 1+jAster}, vecSP, vecLabelsSP );

            end
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            % Order of features : 
            %   Spindle Microtubule (complete optimization)
            %   Astral Microtubules (end optimization)

            vecList = {};
            vecLabelsList = {};
            objList = {};

            % Define properties to get
            props.spindle = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            props.spb = {'amplitude', 'sigma'};
            props.mt = {'endPosition', 'amplitude', 'sigma'};

            % Spindle vector
            [ vecList{1} , vecLabelsList{1}] = getVec( obj.featureList{1}, props.spindle );
            objList{1} = obj.featureList{1};
            numFeatures = 1;

            % SPB 1 vector
            [ vecList{2} , vecLabelsList{2}] = getVec( obj.featureList{2}.featureList{1}, props.spb);
            objList{2} = obj.featureList{2}.featureList{1};
            numFeatures = numFeatures+1;

            % SPB 2 vector
            [ vecList{3} , vecLabelsList{3}] = getVec( obj.featureList{3}.featureList{1}, props.spb);
            objList{3} = obj.featureList{3}.featureList{1};
            numFeatures = numFeatures+1;

            % Microtubule vector : Aster 1
            for jmt = 2 : obj.featureList{2}.numFeatures
                [ vecList{ numFeatures+jmt} , vecLabelsList{ numFeatures+jmt}] = getVec( obj.featureList{2}.featureList{1+jmt}, props.mt );
                objList{ numFeatures+jmt} = obj.featureList{2}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{2}.numFeatures-1;

            % Microtubule vector : Aster 2
            for jmt = 2 : obj.featureList{3}.numFeatures
                [ vecList{ numFeatures+jmt} , vecLabelsList{ numFeatures+jmt}] = getVec( obj.featureList{3}.featureList{1+jmt}, props.mt );
                objList{ numFeatures+jmt} = obj.featureList{3}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{3}.numFeatures-1;

            % Prepend environmental parameters to each vec and vecLabel
            [vecE, vecLabelsE] = getVecEnvironment( obj, obj.props2Fit);
            
            for jFeat = 1 : numFeatures
                vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
                vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
            end

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
        function obj = addSubFeatures( obj, Image2Find)
            % Searches for missing features and adds them to the featureList
            
            % For a spindle, we can only add more astral microtubules
            % For each pole we will find a possible microtubule to add, then we will pick the best of the two possibilities and add them to our feature list.

            % How to find a possible microtubule:
            %   1. Angular Sweep (not clear how to ensure a line is found if there are no peaks)
            %   2. Pick the location of the highest residual and draw a line from both pole connecting to the max residual point. Find the mean intensity of each line. Subtract the mean from the actual voxel values along the line to find a total residual. Pick the pole with the minimum residual.

            props2Fit = {'endPosition', 'amplitude', 'sigma'};
            display = {'Color', 'Red', 'LineWidth', 3};
            if obj.dim==3, sigma=[1.2 1.2 1.0]; elseif obj.dim==2, sigma=[1.2 1.2]; end

            % Get Spindle Angle
            spindleAngle = getAngleXY( obj);
            spindleExclusionRange = deg2rad( 10);

            % ask subfeatures to find missing features
            feature{1} = obj.featureList{2}.findBestMissingFeature( Image2Find, spindleAngle(1), spindleExclusionRange);
            feature{2} = obj.featureList{3}.findBestMissingFeature( Image2Find, spindleAngle(2), spindleExclusionRange);

            % Make higher level decision on the best missing feature
            [~,idxKeep] = min( [ feature{1}.residual, feature{2}.residual ] );

            % Create feature object 
            amp = feature{ idxKeep}.amplitude;
            feature = Line( feature{idxKeep}.startPosition, feature{idxKeep}.endPosition, amp, sigma, obj.dim, props2Fit, display);
%             feature.fillParams();
            
            % Add feature to the correct subfeature 
            idxAdd = obj.featureList{1+idxKeep}.addFeatureToList( feature );
            % Add feature ID and update the featureMap
            feature.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
            obj.featureMap( feature.ID) = [ 1+idxKeep idxAdd];
                
        end
        % }}}

        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % For a spindle, we can only remove  astral microtubules
            % For each pole we will find a possible microtubule to remove, then we will pick the best of the two possibilities and remove them from our feature list.
            % Look at the residual under each astral microtubule
            % Remove the one with the biggest mean residual.
            % If there are no microtubules to remove
            successRemove = 1;

            % Ask asters to give their worst microtubule features
            worstFeature{1} = obj.featureList{2}.findWorstFeature( Image2Find);
            worstFeature{2} = obj.featureList{3}.findWorstFeature( Image2Find);

            % If no features to remove, exit the function
            if worstFeature{1}.idx == 0 && worstFeature{2}.idx == 0
                successRemove = 0;
                return
            end

            % Make higher level decision on the worst-est feature
            [~, idxAster ] = max( [ worstFeature{1}.residual, worstFeature{2}.residual] );
            idxMT = 1 + worstFeature{ idxAster}.idx;

            % Remove the worst microtubule
            obj.featureList{ 1 + idxAster}.removeFeatureFromList( idxMT);
                
        end
        % }}}

        % findEnvironmentalConditions {{{
        function obj = findEnvironmentalConditions( obj)
            % Finds the cytoplasmic and the nucleoplasmic backgrounds
            
            imVals = obj.image( obj.image(:) > 0);
            obj.background = median( imVals);

            % To find the nuclear background, we will raise a threshold value until only

        end
        % }}}

        % getVecEnvironment {{{
        function [vec, vecLabels] = getVecEnvironment( obj, props2get)

            if nargin < 2
                props2get = {'background', 'backgroundNuclear'};
            end

            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
            end

            % Get vector of Properties
            % Also get a string array with property names
            vec = []; vecLabels = {};
            for jProp = 1 : length( props2get)
                vec = [ vec, obj.( props2get{jProp} ) ];
                if numel( obj.( props2get{jProp} ) ) ~= length( obj.( props2get{jProp} ) )
                    error('getVec : property selected is a non-singleton matrix')
                end
                numRep = length( obj.( props2get{jProp} ) );
                labelRep = cell( 1, numRep);
                labelRep(:) = props2get(jProp);
                vecLabels = { vecLabels{:}, labelRep{:} };
            end

        end
        % }}}
        
        % absorbVecEnvironment {{{
        function obj = absorbVecEnvironment( obj, vec, vecLabels)
        
            props2find = {'background', 'background_nuclear'};

            % find the index of start positions
            for jProp = 1 : length( props2find)
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );
                
                % Checking
                if length( idxProp) == 0
                    continue
                end
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVec: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}
        
        % simulateAll {{{
        function imageOut = simulateAll( obj, imageIn, featureID)
           
            imageOut = imageIn;
            % Only simulate feature with matching ID
            % Find the feature
            cFeature = obj.findObjectFromID( featureID);
            imageOut = cFeature.simulateFeature( size(imageIn) );

            % Fill background where features are not prominent-er
            if ~isempty( obj.background)
%                 imageOut( imageOut < obj.background) = obj.background;
                imageOut = imageOut + obj.background;
            end

            % Add nuclear background
            if ~isempty( obj.backgroundNuclear)
                imageOut = imageOut + obj.backgroundNuclear.*obj.maskNuclear;
            end

            % Add Cellular Mask
            imageOut = imageOut .* logical( obj.image);
            
        end
        % }}}
        
        % syncFeaturesWithMap {{{
        function obj = syncFeaturesWithMap( obj)
            % update the feature ids to reflect any initiliazation/additions/removals
            % update the map with the features ids and their locations
            
            map = obj.featureMap;
            counter = uint32(max( cell2mat( keys( map) ) ) + 1);
            if isempty( counter)
                % initialize it
                counter = uint32(1);
            end
            
            % assign unique ID to features
            
            % Depth 0 (self)
            loc = [];
            obj.ID = counter; 
            map( counter) = loc; counter = counter+1;
            
            % Depth 1 (features)
            for jFeat = 1 : obj.numFeatures
                obj.featureList{jFeat}.ID = counter; 
                map( counter) = [ loc jFeat]; counter=counter+1;
            end
            
            % Depth 2 (subfeatures)
            % assign unique ID to subfeatures
            for j1 = 1 : obj.numFeatures
                % if features are organizers
                if isprop( obj.featureList{j1}, 'numFeatures')
                    %assign unique ID to subfeatures
                    for j2 = 1 : obj.featureList{j1}.numFeatures
                        obj.featureList{j1}.featureList{j2}.ID = counter;
                        map( counter) = [ loc j1 j2]; counter=counter+1;
                    end
                end
            end
            
            % Depth 3 (subsubfeatures - overkill)
            % assign unique ID to subsubfeatures
            for j1 = 1 : obj.numFeatures
                % if features are organizers
                if isprop( obj.featureList{j1}, 'numFeatures')
                    for j2 = 1 : obj.featureList{j1}.numFeatures
                        % if subfeatures are organizers
                        if isprop( obj.featureList{j1}.featureList{j2}, 'numFeatures')
                            %assign unique ID to subsubfeatures
                            for j3 = 1 : obj.featureList{j1}.numFeatures
                                obj.featureList{j1}.featureList{j2}.featureList{j3}.ID = counter;
                                map( counter) = [ loc j1 j2 j3]; counter=counter+1;
                            end
                        end                        
                    end
                end
            end
            
        end
        % }}}
        
        % findObjectFromID {{{
        function subObj = findObjectFromID( obj, searchID)
            % use the featureMap to locate and reutrn the object with ID
            % searchID
            
            try
                objLoc = obj.featureMap( searchID);
            catch
                error('findObjectFromID : searchID did not match any ID in featureMap')
            end
            
            switch length( objLoc)
                case 0
                    subObj = obj;
                case 1
                    subObj = obj.featureList{ objLoc(1)};
                case 2
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) };
                case 3
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) }.featureList{ objLoc(3) };
                case 4
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) }.featureList{ objLoc(3) }.featureList{ objLoc(4) };
                otherwise
                    error('findObjectFromID : your object is too deep for this function')
            end
                
        end
        % }}}       

        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.image = uint16( obj.image);
            S.props2Fit = obj.props2Fit;

            for jFeat = 1 : length( obj.featureList)
                S.featureList{ jFeat} = saveAsStruct( obj.featureList{jFeat} );
            end

        end
        % }}}
        
    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Spindle')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = Line.loadFromStruct( S.featureList{ 1} ); 
            for jFeat = 2 : length( S.featureList)
                featureList{ jFeat} = AsterMT.loadFromStruct( S.featureList{ jFeat} ); 
            end

            obj = Spindle( S.dim, S.image, featureList, S.props2Fit);
%             obj = obj@Organizer( S.dim, featureList, S.type);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap()

        end
        % }}}

    end

end
