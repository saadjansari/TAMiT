classdef SpotBank < OrganizerMaster 
% A SpotBank is a higher level feature composed of some number of Spot objects
    properties
    end

    methods (Access = public)

        % SpotBank {{{
        function obj = SpotBank( dim, image, featureList, props2Fit)

            for jspot = 1 : length(featureList)
                if ~strcmp( featureList{jspot}.type, 'Spot')
                    error('SpotBank: spots must be a cell array containing ''Spot'' objects') 
                end
            end
            obj = obj@OrganizerMaster( dim, image, featureList, props2Fit, 'SpotBank');

        end
        % }}}

        % getVec {{{
        function [vec, vecLabels,ub,lb] = getVec( obj, propsSpot)

            if nargin < 2
                propsSpot = {'position', 'amplitude', 'sigma'};
            end
            
            % Get general environmental properties
            [vecE, vecLabelsE] = obj.getVecEnvironment();
            ubE = vecE; lbE = vecE;
            
            vec = []; lb = []; ub = [];
            vecLabels = {};
            % Loop over spots and get their vectors. 
            for js = 1 : obj.numFeatures
                [vec_spot, vecLabels_spot, ub_spot, lb_spot] = getVec( obj.featureList{ js}, propsSpot );
                vec = [vec , vec_spot]; lb = [lb, lb_spot]; ub = [ub, ub_spot];
                vecLabels_spot = strcat( 'Spot', num2str( js), '_', vecLabels_spot);
                vecLabels = { vecLabels{:}, vecLabels_spot{:} };
            end
            if isempty(lb)
                clearvars lb ub
            end
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList, ubList, lbList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects

            vecList = {};
            vecLabelsList = {};
            objList = {};
            ubList = {}; lbList = {};

            % Define properties to get
            propsSpot = {'position', 'amplitude', 'sigma'};

            % Spot vectors :
            for jmt = 1 : obj.numFeatures
                [ vecList{ jmt} , vecLabelsList{ jmt}, ubList{jmt}, lbList{jmt}] = getVec( obj.featureList{jmt}, propsSpot );
                objList{ jmt} = obj.featureList{jmt};
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

            % Get the vector for each spot, and ask the spot objects to absorb the vector
            for js = 1 : obj.numFeatures
                
                strSpot = ['Spot', num2str(js), '_'];
                idxSpot = find( ~cellfun( @isempty, strfind( vecLabels, strSpot ) ) );
                vecSpot = vec( idxSpot );
                vecLabelsSpot = erase( vecLabels( idxSpot), strSpot );
                obj.featureList{ js} = absorbVec( obj.featureList{ js}, vecSpot, vecLabelsSpot, errBoolean);

            end
            
        end
        % }}}

        % addSubFeatures {{{
        function obj = addSubFeatures( obj, image2Find)
            % Searches for missing features and adds them to the featureList
            
            % Kinetochore Bank:
            %   Take the residual image (Image2Find) and find the brightest pixel, this is a new kinetochore
           
            props2Fit = {'position', 'amplitude', 'sigma'};
            
            props = Cell.GetFeatureProps();
            display = props.spot.graphics.red;
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
                
        function obj = finalizeAddedFeatures(obj)
            
%             % Display properties for lines in asters
%             display = {'Color', [0 0.8 0], 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 1};
%             
%             for jL = 1 : length(obj.featureList )
%                 obj.featureList{jL}.display = display;
%             end
                
        end
        
        % }}}
        
    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'SpotBank')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = Spot.loadFromStruct( S.featureList{ 1} ); 
            for jFeat = 1: length( S.featureList)
                featureList{ jFeat} = Spot.loadFromStruct( S.featureList{jFeat} );
            end

            obj = SpotBank( S.dim, S.image, featureList, S.props2Fit);
            obj.findEnvironmentalConditions();
            obj.syncFeaturesWithMap();
        end
        % }}}

    end

end
