classdef KinetochoreBank < Organizer 
% A KinetochoreBank is a higher level feature that sits inside a mitotic cell. It is composed of some number of Kinetochore objects
    properties
        background
        backgroundNuclear
        maskNuclear
        featureMap
        image
    end

    methods (Access = public)

        % KinetochoreBank {{{
        function obj = KinetochoreBank( image, varargin)

            if nargin == 0 
                error( 'KinetochoreBank: number of inputs is inconsistent')
            elseif nargin == 1
                warning('KinetochoreBank: initialized with 0 kinetochores')
            end

            for jKC = 1 : length(varargin)
                if ~strcmp( varargin{jKC}.type, 'Spot')
                    error('KinetochoreBank: featureList must be a cell array containing ''Spot'' objects') 
                end
            end

            dim = length( size( image) );
            obj = obj@Organizer( dim, varargin, 'KinetochoreBank');
            obj.image = image;

            % initialize featureMap and assign IDs
            obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            % assign voxels to its basic elements
            for jFeat = 1 : length( obj.featureList)
                obj.featureList{jFeat}.findVoxelsInsideMask( logical( image) );
            end

        end
        % }}}

        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props)

            vec = [];
            vecLabels = {};

            if nargin < 2
                props.kc = {'position', 'amplitude', 'sigma'};
            end

            % Loop over kinetochores and get their vectors. 
            for jkc = 1 : length( obj.featureList)
                [vec_kc, vecLabels_kc] = getVec( obj.featureList{ jkc}, props.kc );
                vec = [vec , vec_kc];
                vecLabels_kc = strcat( 'KC', num2str( jkc), '_', vecLabels_kc);
                vecLabels = { vecLabels{:}, vecLabels_kc{:} };
            end

        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj, props)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            vecList = {};
            vecLabelsList = {};
            objList = {};

            % Define properties to get
            if nargin < 2
                props.kc = {'position', 'amplitude', 'sigma'};
            end

            % Kinetochore vectors 
            for jkc = 1 : length( obj.featureList )
                [ vecList{ jkc} , vecLabelsList{ jkc}] = getVec( obj.featureList{jkc}, props.kc );
                objList{ jkc} = obj.featureList{jkc};
            end
            numFeatures =  length( obj.featureList);

        end
        % }}}
        
        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)

            % Get the vector for each kinetochore, and ask the kinetochore objects to absorb the vector
            for jkc = 1 : length( obj.featureList)
                
                strKC = ['KC', num2str(jkc), '_'];
                idxKC = find( ~cellfun( @isempty, strfind( vecLabels, strKC ) ) );
                vecKC = vec( idxKC );
                vecLabelsKC = erase( vecLabels( idxKC), strKC );
                obj.featureList{ jkc} = absorbVec( obj.featureList{ jkc}, vecKC, vecLabelsKC);

            end
        end
        % }}}

        % addSubFeatures {{{
        function obj = addSubFeatures( obj, image2Find)
            % Searches for missing features and adds them to the featureList
            
            % Kinetochore Bank:
            %   Take the residual image (Image2Find) and find the brightest pixel, this is a new kinetochore
           
            props2Fit = {'position', 'amplitude', 'sigma'};
            display = {'Color', [0 0.8 0] , 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 1};
            if obj.dim==3, sigma=[1.2 1.2 1.0]; elseif obj.dim==2, sigma=[1.2 1.2]; end

            % Convolve image, then look for max value inside nuclear region
            imageG = imgaussfilt( image2Find, 1) .* obj.maskNuclear;
            image2D = max( imageG, [], 3);

            [amplitude, idx] = max( image2D(:) );
            [y, x] = ind2sub( size(image2D), idx);
            [~, idxMax] = max( image2Find, [], 3);
            z = idxMax( y, x);

            % Create feature object 
            amplitude = image2Find( y,x,z);
            kc = Spot( [ x, y, z], amplitude, sigma, obj.dim, props2Fit, display);
            
            % Add feature to featurelist
            obj.addFeatureToList( kc );
                
        end
        % }}}

        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % Look at the residual under each kinetochore
            % Remove the one with the biggest summed residual.
            successRemove = 1;

            % If no features to remove, exit the function
            if length( obj.featureList) == 0
                successRemove = 0;
                return
            end

            % find the residuals by simulation
            residuals = [];
            for jfeat = 1 : length( obj.featureList)
                % simulate the feature
                imSim = obj.featureList{jfeat}.simulateFeature( size(Image2Find) );
                residuals( jfeat) = sum( abs(Image2Find(:) - imSim(:)).^2);
            end

            % Find the max residual, identify the worst feature and remove it
            [~, idxRm] = max( residuals);
            obj.removeFeatureFromList( idxRm);
                
        end
        % }}}
        
        % findEnvironmentalConditions {{{
        function obj = findEnvironmentalConditions( obj)
            % Finds the cytoplasmic and the nucleoplasmic backgrounds
            
            imVals = obj.image( obj.image(:) > 0);
            obj.background = median( imVals);

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
                searchID
                obj.featureMap.keys
                obj.featureMap.values
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
            S.image = uint16( obj.image);

            for jFeat = 1 : length( obj.featureList)
                S.featureList{ jFeat} = saveAsStruct( obj.featureList{jFeat} );
            end

        end
        % }}}
        
    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'KinetochoreBank')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            for jFeat = 1 : length( S.featureList)
                featureList{ jFeat} = Spot.loadFromStruct( S.featureList{ jFeat} ); 
            end

            obj = KinetochoreBank( S.image, featureList);
%             obj = obj@Organizer( S.dim, featureList, S.type);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap()

        end
        % }}}

    end

end
