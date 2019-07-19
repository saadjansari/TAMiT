classdef BasicElement < Feature

    properties
        amplitude
        sigma
        props2Fit
        display
        params
    end

    methods ( Access = public )

        % BasicElement {{{
        function obj = BasicElement( dim, image, amplitude, sigma, props2Fit, display, type)
        % Basic Element : this is the constructor function for a Basic Element 

            % Ensure dim matches image dimensionality
            if dim ~= length( size( image) )
                error( 'Feature: input argument dim does not match dimensionality of input argument image')
            end

            if nargin < 7
                type = 'generic basic element';
            end

            obj = obj@Feature( dim, image, type);
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
                try
                if ~any( strcmp( props2get{ jProp}, properties( obj) ) )
                    error( 'getVec : unknown property in props')
                end
                catch
                    props2get
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

    end

    methods ( Abstract = true )
       % All subclasses must contain the implemention for the following functions 
        
        % simulateFeature: ability to create a simulated image of a feature
        [ obj, imageSim] = simulateFeature( obj, imageIn)
        
        % displayFeature: ability to display graphically the feature
        ax = displayFeature( obj, ax)

    end

end
