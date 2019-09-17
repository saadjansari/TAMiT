classdef IMTBank < Organizer
% A IMTBank is a higher level feature that sits inside an interphase cell. It is composed of Microtubule asters 
    properties
        numAsters 
        props2Fit
        background
        featureMap
        image
    end

    methods (Access = public)
        
        % IMTBank {{{
        function obj = IMTBank( dim, image, featureList, props2Fit)

            obj = obj@Organizer( dim, featureList, 'IMTBank');
            obj.image = image;
            obj.props2Fit = props2Fit;
            
            % initialize featureMap and assign IDs
            obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            obj.numAsters = length( featureList);

            % assign voxels to its basic elements
            for jAster = 1 : obj.numAsters
                for jFeat = 1 : length( obj.featureList{ jAster}.featureList)

                    obj.featureList{ jAster}.featureList{jFeat}.findVoxelsInsideMask( logical( image) );

                end
            end

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props)

            % Get general environmental properties
            [vec, vecLabels] = getVecEnvironment( obj, obj.props2Fit);

            % get vectors from features
            props.mtoc = {'position', 'amplitude', 'sigma'};
            props.imt = {'endPolyCoef', 'amplitude', 'sigma'};

            % Loop over mt arrays and get their vectors. Remove the MTOCs from the fit vectors. 
            vec = [];
            vecLabels = [];
            for jAster = 1 : obj.numAsters
                [vec_mtarray, vecLabels_mtarray] = getVec( obj.featureList{ jAster}, props.mtoc, props.imt);
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
            for jAster = 1 : obj.numAsters

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

            % Define properties to get
            props.mtoc= {'position', 'amplitude', 'sigma'};
            props.imt = {'endPolyCoef', 'amplitude', 'sigma'};

            % Aster Vector 
            numFeat = 0
            cFeat = 1+numFeat;
            for jAster = 1 : obj.numAsters

                % get mtoc vector
                [ vecList{cFeat} , vecLabelsList{cFeat}] = getVec( obj.featureList{jAster}.featureList{1}, props.mtoc);
                objList{cFeat} = obj.featureList{jAster}.featureList{1};

                % get imt vectors
                nummt = length( obj.featureList{jAster}.featureList) - 1;
                for jmt = 1 : nummt
                    [ vecList{ cFeat+jmt} , vecLabelsList{ cFeat+jmt}] = getVec( obj.featureList{jAster}.featureList{1+jmt}, props.imt );
                    objList{ cFeat+jmt} = obj.featureList{jAster}.featureList{1+jmt};
                end
                numFeat = length( vecList);
                cFeat = 1 : numFeat;
            end

            % Prepend environmental parameters to each vec and vecLabel
            [vecE, vecLabelsE] = getVecEnvironment( obj, obj.props2Fit);
            
            for jFeat = 1 : numFeat
                vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
                vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
            end

        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures( obj)

            % update number of asters 
            obj.numAsters = length( obj.featureList); 

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
            obj.syncFeaturesWithMap();

        end
        % }}}

    end

end
