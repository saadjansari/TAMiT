classdef MT_linear < Feature

    properties
        length
        orientation
        startPosition
        endPosition
    end

    methods
       
        % MT_linear {{{
        function obj = MT_linear( position, amplitude, sigma, dim, ref_image)

            obj = obj@Feature( position, amplitude, sigma, dim, ref_image, 'MT_linear');

            obj.startPosition = obj.position( 1, :);
            obj.endPosition = obj.position( 2, :);

        end
        % }}}

        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props2get)

            % sample props2get
            if nargin==1
            props2get = {'startPosition', 'endPosition', 'amplitude', 'sigma'};
            end
            
            % make sure props input matches properties defined in class
            for jProp = 1 : length( props2get)
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
            end

            % Get vector of Properties
            vec = [];
            for jProp = 1 : length( props2get)
                vec = [ vec, obj.( props2get{jProp} ) ];
            end

            % Also get a string array with property names
            vecLabels = {};
            for jProp = 1 : length( props2get)
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
        
        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)

            props2find = {'startPosition', 'endPosition', 'amplitude', 'sigma'};

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
                obj.( props2find{jProp} );
                vec( idxProp);
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}

        % simulateFeature {{{
        function imageOut = simulateFeature( obj, imageIn)

            if nargin < 2
                imageIn = obj.ref_image;
            end

            % Simulate a gaussian line 
            startPosition = obj.startPosition;
            endPosition = obj.endPosition;
            amplitude = obj.amplitude;
            sigma = obj.sigma;
            imageFeat = 0*imageIn; 
            imageFeat = amplitude*mat2gray( Cell.drawGaussianLine3D( startPosition, endPosition, sigma, imageFeat) );
            imageOut = imageIn;
            imageOut( imageFeat > imageIn) = imageFeat( imageFeat > imageIn);
%             imageOut = imageFeat + imageIn;

        end
        % }}}

        % displayFeature {{{
        function ax = displayFeature( obj, ax)

            if nargin < 2
                f = figure;
                ax = axes; axis ij; hold on;
                imagesc( max( obj.ref_image, [], 3) ); colormap gray; axis equal;
            end

            line( [obj.startPosition(1) obj.endPosition(1)], [obj.startPosition(2) obj.endPosition(2)], 'Color', 'r', 'LineWidth', 2)

        end
        % }}}

    end
end
