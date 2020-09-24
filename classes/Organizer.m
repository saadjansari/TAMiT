classdef Organizer < Feature
% This is an organizer. This is a collection of features that are somehow connected together either spatially or through the fact that share parameters, or both. 
    properties
        featureList
        numFeatures
        featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any')
    end

    methods

        % Organizer {{{
        function obj = Organizer( dim, featureList, type)
        % Organizer : this is the constructor function for an Organizer

            if nargin < 3 
                type = 'generic organizer';
            end

            obj = obj@Feature( dim, type);
            obj.featureList = featureList;
            obj.numFeatures = length( featureList);

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, varargin)

            if nargin < 2
                error('displayFeature: pass in an axes handle to display Feature')
            end

            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                ax = obj.featureList{jFeat}.displayFeature( ax, varargin{:});
            end

        end
        % }}}
        
        % displayFeature3D {{{
        function ax = displayFeature3D( obj, ax,sizeZ)

            if nargin < 2
                error('displayFeature: pass in an axes handle to display Feature')
            end

            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                ax = obj.featureList{jFeat}.displayFeature3D( ax,sizeZ);
            end

        end
        % }}}
        
        % displayFeature {{{
        function ax = displayFeatureXZ( obj, ax)

            if nargin < 2
                error('displayFeature: pass in an axes handle to display Feature')
            end

            if obj.dim == 2
                error('displayFeatureXZ: organizer must be 3-dimensional')
            end
            
            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                ax = obj.featureList{jFeat}.displayFeatureXZ( ax);
            end

        end
        % }}}
        
        % displayFeature {{{
        function ax = displayFeatureYZ( obj, ax)

            if nargin < 2
                error('displayFeature: pass in an axes handle to display Feature')
            end
            if obj.dim == 2
                error('displayFeatureYZ: organizer must be 3-dimensional')
            end
            
            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                ax = obj.featureList{jFeat}.displayFeatureYZ( ax);
            end

        end
        % }}}
        
        % makeCopyDeep {{{
        function objCopy = copyDeep( obj)
            % makes a deep copy of the handle object

            % first make a shallow copy
            objCopy = copyShallow( obj);
            
            % copy featureMap if present
            if isprop( objCopy, 'featureMap') && ~isempty( obj.featureMap.keys )
                objCopy.featureMap = containers.Map( obj.featureMap.keys,obj.featureMap.values);
            else
                obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
            end

            % now make it deep by copying obj.featureList
            for jFeat = 1 : obj.numFeatures

                % if feature is an organizer itself, then do an iterative deepcopy
                if isprop( objCopy.featureList{jFeat}, 'numFeatures' )
                    objCopy.featureList{jFeat} = copyDeep( obj.featureList{jFeat} );
                else
                    objCopy.featureList{jFeat} = copyShallow( obj.featureList{jFeat} );
                end

            end

        end
        % }}}

        % simulateFeature {{{ 
        function [imageOut,error_code, err] = simulateFeature( obj, sizeImage)

            if nargin < 2
                error('simulateFeature: must pass an imageSize to make a simulated Image')
            end
            error_code = 0; err=0;
            % Simulate all the features
            imageOut = zeros( sizeImage);
            for jFeat = 1 : obj.numFeatures
                try
                    [imNew, error_code_new, errr] = simulateFeature( obj.featureList{ jFeat}, sizeImage);
                catch
                    imNew = simulateFeature( obj.featureList{ jFeat}, sizeImage);
                    error_code_new = 0;
                    errr = 0;
                end
                imageOut = imageOut + imNew;
                error_code = any( [error_code, error_code_new]);
                err = max( [err, errr]);
            end
            obj.imageSim = imageOut;

        end
        % }}}

        % addFeatureToList {{{
        function idx = addFeatureToList( obj, featNew)

            % Check to ensure MT is an object of type microtubule
%             if ~strcmp( MT.type, 'MT_linear')
%                 error('addFeature: feature object to add must be of the correct type'); end

            % Add MT
            obj.featureList = { obj.featureList{:}, featNew};
            obj.numFeatures = obj.numFeatures + 1;
            idx = obj.numFeatures;
            obj.syncFeaturesWithMap();

        end
        % }}}

        % removeFeatureFromList {{{
        function obj = removeFeatureFromList( obj, idxFeature) 

            % Get ID of feature
            id = cellfun( @(x) x.ID, obj.featureList(idxFeature));
            
            % Remove feature and its record from the map
            obj.featureList( idxFeature) = [];
            remove(obj.featureMap, id);
            
            % Remove
            obj.numFeatures = obj.numFeatures - length(id);
%             obj.syncFeatures();

        end
        % }}}
        
        % getSubFeatureNumber {{{
        function numFeat = getSubFeatureNumber( obj)
            
            numFeat = 0;
            for jFeat = 1 : obj.numFeatures
                if isprop( obj.featureList{jFeat}, 'featureList')
                    numFeat = numFeat + obj.featureList{jFeat}.getSubFeatureNumber();
                else
                    numFeat = numFeat + 1;
                end
            end

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
                subObj = [];
                warning('findObjectFromID : searchID did not match any ID in featureMap')
                return
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

        % fillparams {{{
        function obj = fillParams( obj, sizeImage)
           
            % fill params for each object
            for jobj = 1 : obj.numFeatures
                obj.featureList{jobj}.fillParams(sizeImage);
            end
            
        end
        % }}}

        % updateFeatureIDs {{{
        function obj = updateFeatureIDs( obj)
            % update the feature ids to reflect any initiliazation/additions/removals
            
            global COUNTER             
            % Assign IDs to its subfeatures
            for jFeature = 1: obj.numFeatures
                obj.featureList{ jFeature}.ID = COUNTER;
                COUNTER = COUNTER + 1;
            end

            % Prompt subfeatures to do the same if they are organizers
            for jFeature = 1: obj.numFeatures

                % Determine if subfeature is an organizer
                feat = obj.featureList{jFeature};
                feat_classname = class( feat);
                feat_superclassnames = superclasses( feat_classname); 

                % If organizer, Sync subfeatures
                if strcmp( feat_classname, 'Organizer') || any( strcmp( feat_superclassnames, 'Organizer') )
                    feat.updateFeatureIDs();
                end

            end
            
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
%             map( counter) = loc; 
            counter = counter+1;
            
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
        
        % updateFeatureMap {{{
        function obj = updateFeatureMap( obj)
            % construct a new feature map from IDs of existing features

            % Clear current feature map
            mapLocal = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            % Update featureMap
            for jFeature = 1: obj.numFeatures
                
                % Determine if subfeature is an organizer
                feat = obj.featureList{jFeature};
                feat_classname = class( feat);
                feat_superclassnames = superclasses( feat_classname);

                mapLocal( feat.ID) = [jFeature];
                
                % If organizer, also append sublocal feature map
                if strcmp( feat_classname, 'Organizer') || any( strcmp( feat_superclassnames, 'Organizer') )
                    
                    % update sublocal feature map
                    feat.updateFeatureMap();
                    mapSubLocal = feat.featureMap;
                    keysSubLocal = cell2mat( keys( mapSubLocal) );
                    for key = keysSubLocal
                        mapLocal( key) = [jFeature mapSubLocal( key)];
                    end
                end

            end
            obj.featureMap = mapLocal;
            
        end
        % }}}
        
        % syncFeatures {{{
        function obj = syncFeatures( obj)
            % update the feature ids to reflect any initiliazation/additions/removals
            % update the map with the features ids and their locations
            
            obj.updateFeatureIDs();
            obj.updateFeatureMap();
            
        end
        % }}}

        % GetProjection2D {{{
        function obj = GetProjection2D( obj)
           % Get 2D projection of feature
           
           % Change object dimensionality
           if obj.dim == 2
               warning('object dimensionality is already 2')
           end
           
           obj.dim = 2;
           
           % Specific object projection
           try
               obj = obj.GetProjection2DSpecific();
           end
           
           % Repeat for subfeatures
           for jFeature = 1: obj.numFeatures
               obj.featureList{jFeature}.GetProjection2D();
           end
            
        end
        % }}}
        
        % GetProjection3D {{{
        function obj = GetProjection3D( obj)
           % Get 2D projection of feature
           
           % Change object dimensionality
           if obj.dim == 3
               warning('object dimensionality is already 3')
           end
           
           obj.dim = 3;
           
           % Specific object projection
           try
               obj = obj.GetProjection3DSpecific();
           end
           
           % Repeat for subfeatures
           for jFeature = 1: obj.numFeatures
               obj.featureList{jFeature}.GetProjection3D();
           end
            
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            % Accomplish 3 tasks
            %   1. Update 3D features from 2D features with matching ID
            %   2. Remove 3D features if 2D counterparts with matching ID
            %   are not present
            %   3. Add 3D features if there are extra 2D features present
            
            % 1. For each 3D feature, get its ID, then look for matching 2D
            % feature with the same ID
            try
                obj.Update3DFrom2DSpecific( obj2D);
            end
            obj.fit = 'zonly';
            
            rmFeat = [];
            for jF = 1 : obj.numFeatures
                
                ID = obj.featureList{jF}.ID;
                feat2D = obj2D.findObjectFromID( ID);
                
                % Mark features whose 2D counterparts arent present
                if isempty(feat2D)
                    rmFeat = [rmFeat, jF];
                else
                    obj.featureList{jF}.Update3DFrom2D( feat2D );
                end
            end
            
            % 3. Create 3D features if there are extra 2D features present
            for jF = 1 : obj2D.numFeatures
                
                ID = obj2D.featureList{jF}.ID;
                feat3D = obj.findObjectFromID( ID);
                
                if isempty( feat3D)
                    % Change feature to 3D
                    feat3D = obj2D.featureList{jF}.GetProjection3D();
                    feat3D.fit = 'zonly';
                    obj.addFeatureToList( feat3D );
                end
                
            end
            
            % 2. Remove 3D features if 2D counterparts not present
            if ~isempty(rmFeat)
                obj.removeFeatureFromList( rmFeat);
            end
            
        end
        % }}}
        
        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
            % Force features to lie inside the mask. This will
            % shorten any curves whose mtoc is within the mask.
            
            % Ask features to force their subfeatures
            for jF = 1 : obj.numFeatures
                obj.featureList{jF}.forceInsideMask( mask);
            end
            
        end
        % }}}

        % GetStructBasicFeatures {{{
        function feats = GetStructBasicFeatures(obj)
            
            feats = [];
            for jf = 1 : obj.numFeatures
                type = obj.featureList{jf}.type;
                if any(strcmp( superclasses( type), 'BasicElement'))
                    feats = Methods.combineStructsIntoArray( feats, obj.featureList{jf}.GetStructInfo() );
                elseif any(strcmp( superclasses( type), 'Organizer'))
                    feats = Methods.combineStructsIntoArray( feats, obj.featureList{jf}.GetStructBasicFeatures() );
                end
            end

        end
        % }}}
        
        function obj = preOptimize(obj, imOrg, imBkg)
           
            % optimize basic objects
            for jFeat = 1 : length( obj.featureList)
                try
                    obj.featureList{jFeat}.preOptimize(imOrg, imBkg);
                end
            end
            
        end

        function imfMask = getMaskWithoutFeatures(obj, img)
           % Obtain a mask with the features removed
                           
            imfMask = ones( size(img));

            % Simulate features one at a time
            for jb = 1 : obj.numFeatures
                imf = obj.featureList{jb}.simulateFeature( size(img) );
                imfMask = imfMask .* (imf < 0.2*max(imf(:) ));
            end
                   
        end
    end

end
