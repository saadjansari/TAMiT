classdef BasicElement < Feature

    properties
        amplitude
        sigma
        props2Fit
        display
        params = [];
        bounds
    end

    methods ( Access = public )

        % BasicElement {{{
        function obj = BasicElement( dim, amplitude, sigma, props2Fit, display, type)
        % Basic Element : this is the constructor function for a Basic Element 

            if nargin < 6 
                type = 'generic basic element';
            end

            obj = obj@Feature( dim, type);
            obj.amplitude = amplitude;
            obj.sigma = sigma;
            obj.props2Fit = props2Fit;
            obj.display = display;

        end
        % }}}

        % getVec {{{
        function [vec, vecLabels] = getVec( obj, props2get)

            if nargin < 2
                props2get = obj.props2Fit;
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
        function obj = absorbVec( obj, vec, vecLabels, props2find)

            if nargin < 4
                props2find = obj.props2Fit;
            end

            % find the index of start positions
            for jProp = 1 : length( props2find)
                idxProp = find( strcmp( props2find{ jProp} , vecLabels) );

                % If property is not free to change, then continue to next iteration
                if isempty( idxProp)
                    continue;
                end
                
                % Checking
                if length( obj.( props2find{ jProp} ) ) ~= length( vec(idxProp) )
                    error( 'absorbVec: length of vector props to absorb does not match the old property size')
                end
            
                % Set final property
                obj.( props2find{ jProp} ) = vec( idxProp);

            end

        end
        % }}}

        % findVoxelsInsideMask {{{
        function obj = findVoxelsInsideMask( obj, mask)
            
            obj.params.idx = find( mask);
            [obj.params.y, obj.params.x, obj.params.z] = ind2sub( size( mask), obj.params.idx);

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
           obj.sigma = obj.sigma(1:2);
           
           % Specific object projection
           try
               obj = obj.GetProjection2DSpecific();
           end
            
        end
        % }}}
        
        % GetProjection3D {{{
        function obj = GetProjection3D( obj)
           % Get 3D projection of feature
           
           % Change object dimensionality
           if obj.dim == 3
               warning('object dimensionality is already 3')
           end
           
           obj.dim = 3;
           obj.sigma(3) = mean( obj.sigma);
           
           % Specific object projection
           try
               obj = obj.GetProjection3DSpecific();
           end
            
        end
        % }}}
        
        % Update3DFrom2D {{{
        function obj = Update3DFrom2D(obj, obj2D)
            
        end
        % }}}

        % forceInsideMask {{{
        function obj = forceInsideMask( obj, mask)
        end
        % }}}

    end

    methods ( Abstract = true )
       % All subclasses must contain the implemention for the following functions 
        
        % simulateFeature: ability to create a simulated image of a feature
        [ obj, imageSim] = simulateFeature( obj, imageIn)
        
        % displayFeature: ability to display graphically the feature
        ax = displayFeature( obj, ax)

    end

end
