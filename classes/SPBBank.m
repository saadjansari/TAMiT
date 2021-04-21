classdef SPBBank < SpotBank 
% A SPBBank is a specialization of a SpotBank with the added restriction of
% the number of spots.
    properties
    end

    methods (Access = public)

        % SPBBank {{{
        function obj = SPBBank( dim, image, featureList, props2Fit)
            obj = obj@SpotBank( dim, image, featureList, props2Fit);
        end
        % }}}

        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, image2Find)
            % Searches for missing features and adds them to the featureList
            
            successAdd = 0;
            % Take the residual image (Image2Find) and find missing spots,
            % up to a total number of 2 spots.
            n_total = 2;
            if obj.numFeatures < n_total
                
                % Get properties
                props = Cell.GetFeatureProps();
                
                % Initilaize standard deviation
                if obj.dim==3, sigma=[1.2 1.2 1.0]; elseif obj.dim==2, sigma=[1.2 1.2]; end

                % Find spots
                [ spb, intBkg] = MitoticCell.findSPB_sid1( image2Find);
                % Consider the best spot
                if length(spb) == 0
                    return
                elseif length(spb) > 1
                    [~,idxKeep] = sort( [spb.amplitude], 'descend');
                    spb = spb( idxKeep (1) );
                    
                    SPB =  Spot( [ spb.x, spb.y, spb.z], spb.amplitude, ...
                        sigma, obj.dim, props.spot.fit{obj.dim}, props.spot.graphics.red);

                    % Add feature to featurelist
                    imbkg = 0*obj.image; imbkg(:) = obj.background;
                    SPB.preOptimize(obj.image, imbkg);
                    obj.addFeatureToList( SPB );
                    
                    % Add feature ID and update the featureMap
                    SPB.ID = max( cell2mat( keys( obj.featureMap) ) )+1;
                    obj.featureMap( SPB.ID) = [ length(obj.featureList) ];
                    successAdd = 1;
                end
                
            end
                
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
        
    end

end
