classdef KinetochoreBank < Feature
% A KinetochoreBank is a higher level feature that sits inside a mitotic cell. It is composed of some number of Kinetochore objects
    properties
    end

    methods (Access = public)

        % KinetochoreBank {{{
        function obj = KinetochoreBank( varargin)

           if nargin == 0 
                error( 'KinetochoreBank: number of inputs is inconsistent')
            end
            
            dim = length( size( varargin{1}.ref_image) );
            obj = obj@Feature( [], [], [], dim, varargin{1}.ref_image, 'KinetochoreBank');
            obj.featureList = { varargin{:} };
            obj.numKinetochores = length( varargin );

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
        function [ vecList, vecLabelsList, objList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            vecList = {};
            vecLabelsList = {};
            objList = {};

            % Define properties to get
            props.kc= {'position', 'amplitude', 'sigma'};

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
            
            [ image2D, idxMax] = max( image2Find, [], 3);
%             image2D = imgaussfilt( image2D, 1):

            [amplitude, idx] = max( image2D(:) );
            [y, x] = ind2sub( size(image2D), idx);
            z = idxMax( y, x);

            % Create feature object 
            dim = length( size( image2Find) ); 
            if dim==2, sigma=[1.2 1.2]; elseif dim==3, sigma=[1.2 1.2 1.0]; end
            kc = Kinetochore( [ x, y, z], amplitude, sigma, dim, obj.ref_image);
            
            % Add feature to featurelist
            obj.addFeature( kc );
                
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
                imSim = obj.featureList{jfeat}.simulateFeature();
                residuals( jfeat) = sum( abs(Image2Find(:) - imSim(:)).^2);
            end

            % Find the max residual, identify the worst feature and remove it
            [~, idxRm] = max( residuals);
            obj.removeFeature( idxRm);
                
        end
        % }}}
        
    end

end
