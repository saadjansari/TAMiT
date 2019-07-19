classdef Organizer < Feature
% This is an organizer. This is a collection of features that are somehow connected together either spatially or through the fact that share parameters, or both. 
    properties
        featureList
        numFeatures
    end

    methods

        % Organizer {{{
        function obj = Organizer( dim, image, featureList, type)
        % Organizer : this is the constructor function for an Organizer
            
            % Ensure dim matches image dimensionality
            if dim ~= length( size( image) )
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end

            if nargin < 4
                type = 'generic organizer';
            end

            obj = obj@Feature( dim, image, type);
            obj.featureList = featureList;
            obj.numFeatures = length( featureList);

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.image, [], 3) ); colormap gray; axis equal;
            end

            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                ax = obj.featureList{jFeat}.displayFeature( ax);
            end

        end
        % }}}
        
        % makeCopyDeep {{{
        function objCopy = copyDeep( obj)
            % makes a deep copy of the handle object

            % first make a shallow copy
            objCopy = copyShallow( obj);
            
            % copy featureMap if present
            if isprop( objCopy, 'featureMap')
                objCopy.featureMap = containers.Map( obj.featureMap.keys,obj.featureMap.values);
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
        function imageOut = simulateFeature( obj, imageIn)

            if nargin < 2
                imageIn = 0*obj.image;
            end
            
            imageOut = imageIn;

            % Simulate all the features
            for jFeat = 1 : obj.numFeatures
                imageOut = simulateFeature( obj.featureList{ jFeat}, imageOut);
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

        end
        % }}}

        % removeFeatureFromList {{{
        function obj = removeFeatureFromList( obj, idxFeature) 

            obj.featureList( idxFeature) = [];
            obj.numFeatures = obj.numFeatures - 1;

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
        % Looks at the objects inside obj.featureList and returns the matching object            
    
            % Can i get a list of object ID's at this level? For organizers, i should probe them to search for IDs at their level
            subObj = [];

            % Is this object the one?
            if obj.ID == searchID
                subObj = obj;
                return
            end

            % What about its subfeatures?
            for jFeat = 1 : obj.numFeatures
                if obj.featureList{jFeat}.ID == searchID
                    subObj = obj.featureList{jFeat};
                end
            end

            % What about its subsubfeatures?
            if isempty( subObj)
                % check if subfeatures contains subsubfeatures
                for jFeat = 1 : obj.numFeatures
                    if isprop( obj.featureList{jFeat}, 'featureList')
                        subObj = findObjectFromID( obj.featureList{jFeat}, searchID);
                    end
                end
            end
                    
        end
        % }}}

        % fillParams {{{
        function obj = fillParams( obj)
           
            % fill params for each object
            for jObj = 1 : obj.numFeatures
                obj.featureList{jObj}.fillParams();
            end
            
        end
        % }}}

    end

end
