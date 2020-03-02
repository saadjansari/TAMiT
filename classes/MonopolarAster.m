classdef MonopolarAster < Organizer
% A MonopolarAster is a higher level feature that sits inside a monopolar cell. It is composed of 1 AsterMT 
    properties
        props2Fit
        background
        backgroundNuclear 
        maskNuclear
        image
    end

    methods (Access = public)
        
        % MonopolarAster {{{
        function obj = MonopolarAster( dim, image, aster, props2Fit)

            % ensure that there is only 1 element in featureList and that is of type AsterMT
            if length( aster) ~= 1 || ~strcmp( aster{1}.type, 'AsterMT')
                error('MonopolarAster: aster must be of type ''AsterMT'' ') 
            end

            obj = obj@Organizer( dim, aster, 'MonopolarAster');
            obj.image = image;
            obj.props2Fit = props2Fit;
            
            % initialize featureMap and assign IDs
            %obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            % assign voxels to its basic elements
            for jFeat = 1 : length( obj.featureList{1}.featureList)
                obj.featureList{1}.featureList{jFeat}.findVoxelsInsideMask( logical( image) );
            end

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props)

            % Get general environmental properties
            [vec, vecLabels] = getVecEnvironment( obj, obj.props2Fit);

            % get vectors from features
            props.spb = {'position', 'amplitude', 'sigma'};
            props.mt = {'endPosition', 'amplitude', 'sigma'};

            % get vectors from aster. Remove the SPBs from the fit vectors. 
            [vec, vecLabels] = getVec( obj.featureList{1}, props.spb, props.mt);

        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels);

            % Aster Vector
            obj.featureList{ 1}.absorbVec( vec, vecLabels );
            
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
            props.spb = {'position', 'amplitude', 'sigma'};
            props.mt = {'endPosition', 'amplitude', 'sigma'};

            % SPB vector
            [ vecList{1} , vecLabelsList{1}] = getVec( obj.featureList{1}.featureList{1}, props.spb);
            objList{1} = obj.featureList{1}.featureList{1};
            numFeatures = 1;

            % Microtubule vectors : Aster 
            for jmt = 1 : obj.featureList{1}.numFeatures-1
                [ vecList{ numFeatures+jmt} , vecLabelsList{ numFeatures+jmt}] = getVec( obj.featureList{1}.featureList{1+jmt}, props.mt );
                objList{ numFeatures+jmt} = obj.featureList{1}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{1}.numFeatures-1;

            % Prepend environmental parameters to each vec and vecLabel
            [vecE, vecLabelsE] = getVecEnvironment( obj, obj.props2Fit);
            
            for jFeat = 1 : numFeatures
                vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
                vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
            end

        end
        % }}}

        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, Image2Find)
            % Searches for missing features and adds them to the featureList
            
            % For a monopolar aster, we can only add more astral microtubules : we will find a possible microtubule to add 

            % How to find a possible microtubule:
            %   1. Angular Sweep (not clear how to ensure a line is found if there are no peaks)
            %   2. Pick the location of the highest residual and draw a line from both pole connecting to the max residual point. Find the mean intensity of each line. Subtract the mean from the actual voxel values along the line to find a total residual. Pick the pole with the minimum residual.

            props2Fit = {'endPosition', 'amplitude', 'sigma'};
            display = {'Color', [1 0.5 0], 'LineWidth', 3};
            if obj.dim==3, sigma=[1.2 1.2 1.0]; elseif obj.dim==2, sigma=[1.2 1.2]; end

            % ask subfeatures to find missing features
            featureKeep = obj.featureList{1}.findBestMissingFeature( Image2Find, [], []);

            % Create feature object 
            amp = featureKeep.amplitude;
            feature = Line( featureKeep.startPosition, featureKeep.endPosition, amp, sigma, obj.dim, props2Fit, display);
            
            % Add feature to the correct subfeature 
            idxAdd = obj.featureList{1}.addFeatureToList( feature);
            successAdd = 1;
            
        end
        % }}}

        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % For a monopolar aster, we can only remove astral microtubules
            successRemove = 1;

            % Ask asters to give their worst microtubule features
            worstFeature = obj.featureList{1}.findWorstFeature( Image2Find);

            % If no features to remove, exit the function
            if worstFeature.idx == 0 && worstFeature{2}.idx == 0
                successRemove = 0;
                return
            end

            % Remove the worst microtubule
            idxMT = 1+worstFeature.idx;
            obj.featureList{ 1}.removeFeatureFromList( idxMT);
                
        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures(obj)
            
            % ensure aster pole is same as start positions of all microtubules
            aster = obj.featureList{1};
            pole = aster.featureList{1};
            pos = pole.position;

            for jmt = aster.numFeatures-1
                aster.featureList{jmt+1}.startPosition = pos;
            end

        end
        % }}}
        
        function obj = finalizeAddedFeatures(obj)
            
            % Display properties for lines in asters
            display = {'Color', [1 0 1], 'LineWidth', 3};
            
            for jL = 2 : length(obj.featureList{1} )
                obj.featureList{1}.featureList{jL}.display = display;
            end
                
        end

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

            % If any features are outside cellular mask, blow up
            imOut = imageOut .* ~(logical(obj.image));
            if max(imOut(:)) > obj.background + 0.02*max(imageOut(:))
                imageOut(:) = Inf;
            end
            
            % Add Cellular Mask
            imageOut = imageOut .* logical( obj.image);
            
        end
        % }}}
        
        % findObjectFromID {{{
        function subObj = findObjectFromID( obj, searchID)
            % use the featureMap to locate and reutrn the object with ID
            % searchID
            
            % Check Current Object ID
            if searchID == obj.ID
                subObj = obj;
                return
            end
            
            % Check SubFeatures
            try
                objLoc = obj.featureMap( searchID);
            catch
                error('findObjectFromID : searchID did not match any ID in featureMap')
            end
            
            switch length( objLoc)
                case 1
                    subObj = obj.featureList{ objLoc(1)};
                case 2
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) };
                case 3
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) }.featureList{ objLoc(3) };
                case 4
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) }.featureList{ objLoc(3) }.featureList{ objLoc(4) };
                case 5
                    subObj = obj.featureList{ objLoc(1)}.featureList{ objLoc(2) }.featureList{ objLoc(3) }.featureList{ objLoc(4) }.featureList{ objLoc(5) };
                otherwise
                    error('findObjectFromID : your object is too deep for this function. And jesus loves you.')
            end
                
        end
        % }}}       

        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.image = im2uint16( obj.image);
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
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'MonopolarAster')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = AsterMT.loadFromStruct( S.featureList{ 1} ); 

            obj = MonopolarAster( S.dim, S.image, featureList, S.props2Fit);
%             obj = obj@Organizer( S.dim, featureList, S.type);
            obj.findEnvironmentalConditions();
            %obj.syncFeaturesWithMap();

        end
        % }}}

    end

end
