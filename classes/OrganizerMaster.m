classdef OrganizerMaster < Organizer
    % This is a master organizer. There will be implementations of this in the highest level organizers (containing images).
    properties
        image
        mask
        background
        backgroundNuclear
        maskNuclear
        props2Fit
    end

    methods( Access = public)

        % OrganizerMaster {{{
        function obj = OrganizerMaster( dim, image, featureList, props2Fit, type)
        % OrganizerMaster : this is the constructor function for an OrganizerMaster

            if nargin < 5 
                type = 'generic master organizer';
            end

            obj = obj@Organizer( dim, featureList, type);
            obj.image = image;
            obj.props2Fit = props2Fit;
            obj.mask = logical( obj.image ~= 0);

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax, col3D)

            if nargin < 2
                error('displayFeature: pass in an axes handle to display Feature')
            end
            
            % Ask subfeatures to display themselves
            for jFeat = 1 : obj.numFeatures
                if nargin == 3 && col3D == 1
                    ax = obj.featureList{jFeat}.displayFeature( ax, size(obj.image,3));
                else
                    ax = obj.featureList{jFeat}.displayFeature( ax);
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
                imVals = obj.image( obj.image(:) > 0);
                obj.background = median( imVals);
            end

            % Nuclear Background
            if isBkgNuc 
                warning( [obj.type, ': background nuclear search not set up yet.'])
            end

        end
        % }}}

        % getVecEnvironment {{{
        function [vec, vecLabels] = getVecEnvironment( obj, props2get)

            if nargin < 2
                props2get = obj.props2Fit.fit{obj.dim}.Environment;
            end

            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVecEnvironemnt : unknown property in props')
                end
            end

            % Get vector of Properties
            % Also get a string array with property names
            vec = []; vecLabels = {};
            for jProp = 1 : length( props2get)
                vec = [ vec, obj.( props2get{jProp} ) ];
                if numel( obj.( props2get{jProp} ) ) ~= length( obj.( props2get{jProp} ) )
                    error('getVecEnvironemnt: property selected is a non-singleton matrix')
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
        
            props2find = obj.props2Fit.fit{obj.dim}.Environment;

            % find the index of start positions
            for jProp = 1 : length( props2find)
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );
                
                % Checking
                if length( idxProp) == 0
                    continue
                end
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVecEnvironment: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}

        % simulateAll {{{
        function [imageOut, err] = simulateAll( obj, imageIn, featureID)
           
            if ~strcmp(featureID, 'env')

                % Find the feature to simulate using ID
                cFeature = obj.findObjectFromID( featureID);

                % Simulate feature
                [imageOut,error_code, err1] = cFeature.simulateFeature( size(imageIn) );

            else % Just get image
                if isempty( obj.imageSim)
                    % Simulate feature 
                    try
                        [imageOut,error_code, err1] = obj.simulateFeature( size(imageIn) );
                        obj.imageSim = imageOut;
                    catch
                        imageOut = obj.simulateFeature( size(imageIn) ); err1 = 0; error_code =0;
                    end
                else
                    imageOut = obj.imageSim; err1 = 0; error_code =0;
                end

            end
            
            % Penalize if features exceed mask region
            err2 = PenalizeOutsideMask( imageOut, obj.mask);
            
            %if error_code
                %imageOut = ones(size(imageOut));
                %return
            %end
            err = 1+err1+err2;
            
            % Fill background where features are not prominent-er
            if ~isempty( obj.background)
                imageOut = imageOut + obj.background;
            end

            % Add nuclear background
            if isprop(obj, 'backgroundNuclear') && ~isempty( obj.backgroundNuclear)
                imageOut = imageOut + obj.backgroundNuclear.*obj.maskNuclear;
            end

            % Add Cellular Mask
            imageOut = imageOut .* obj.mask;

            % PenalizeOutsideMask {{{
            function amt = PenalizeOutsideMask( imageIn, mask)
                % Penalize features outside 2D mask

                if nargin < 2
                    mask = obj.mask;
                end
               
                % Check if escape mask
                [outside, amt] = Methods.CheckEscapeMask( imageIn, mask);
%                 scaling = 0.02;
%                 if outside
%                     imageOut = imageIn*(1+scaling*0.1) -amt*scaling*obj.image;
%                     disp('   Warning: OrganizerMaster.simulateAll - feature escaped 2D mask, applying scaled cost to force reflection')
%                 else
%                     imageOut = imageIn;
%                 end
            end
            % }}}
            
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
     
        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.dim = obj.dim;
            S.image = im2double( obj.image);
            S.props2Fit = obj.props2Fit;
            S.mask = obj.mask;

            for jFeat = 1 : length( obj.featureList)
                S.featureList{ jFeat} = saveAsStruct( obj.featureList{jFeat} );
            end

        end
        % }}}

        % GetProjection2DSpecific {{{
        function obj = GetProjection2DSpecific( obj)
           % Get 2D projection of feature
           
           % Check object dimensionality
           if obj.dim == 2
               warning('object dimensionality is already 2')
           end
           
           obj.image = max( obj.image, [], 3);          
           try
               obj.mask = min( obj.mask, [], 3);  
           end
           
        end
        % }}}

        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
            % Force features to lie inside the mask. This will
            % shorten any curves whose mtoc is within the mask.
            % This will remove any asters whose mtoc is outside the mask
            
            if nargin < 2
                mask = obj.mask;
            end
            
            % Ask features to force their subfeatures
            %rmFeat = [];
            for jF = 1 : obj.numFeatures
                obj.featureList{jF}.forceInsideMask( mask);
                %if success
                    %rmFeat = [ rmFeat, jF];
                %end
            end

            % Remove bad features
            %obj.removeFeaturesFromList( rmFeat);
            
        end
        % }}}
        
        % updateSubFeatures {{{
        function obj = updateSubFeatures(obj)
            obj.numFeatures = length( obj.featureList);
        end
        % }}}

        function obj = finalizeAddedFeatures(obj)
        end

        function obj = preOptimize(obj)
           
            % optimize basic objects
            imbkg = 0*obj.image; imbkg(:) = obj.background;
            for jFeat = 1 : length( obj.featureList)
                try
                    obj.featureList{jFeat}.preOptimize( obj.image, imbkg);
                end
            end
            
        end
        
    end
end



