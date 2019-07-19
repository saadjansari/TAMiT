classdef Microtubule < Feature
% This is a specialization of the Feature superclass for a microtubule
    properties
        polyForm % is the microtubule represented by polynomial coefficients (boolean)
        polyOrder % order of the polynomial. this is a 1 x dim vector containing polynomial order in the X, Y (and possibly Z) dimensions
        polyCoef % this is a structure with fields 'X' and 'Y', and possibly 'Z', which contains the polynomial order
        startPosition
        endPosition
    end

    methods
        
        % Microtubule {{{
        function obj = Microtubule( position, amplitude, sigma, dim, ref_image, polyForm) 
        % This is the constructor function for the Microtubule specialization of a Feature
            
            % Call the superclass Feature constructor
            obj = obj@Feature( position, amplitude, sigma, dim, ref_image, 'Microtubule');

            obj.startPosition = obj.position( 1, :);
            obj.endPosition = obj.position( 2, :);

            % default obj.polyForm to 0 if not provided
            % nargin
            if nargin < 6, obj.polyForm = 0;
            else, obj.polyForm = polyForm; end

        end
        % }}}

        function vec = getVec( props2get)



        end

        function simulatedFeature = simulateFeature( obj )
        % simulateFeature: simulates the feature 
        end

        function featureHandle = displayFeature( obj, dimDisplay, displayProps)
        % displayMicrotubule : calls either displayMicrotubule2D or displayMicrotubule3D
        end

        function imageNew = addSimulatedFeature( obj, imageOld)
        end

        % getPropertiesToFit {{{
        function fitProperties = getPropertiesToFit( obj, props)
        % getPropertiesToFit : props is a string array whose strings must case-match the properties of the class.

            % make sure props input matches properties defined in class
            for jProp = 1 : length( props)
                if ~any( strcmp( props{ jProp}, properties( obj) ) )
                    error( 'getPropertiesToFit : unknown property in props')
                end
            end

            % Get structure of Properties
            for jProp = 1 : length( props)
                structProperties.( props{jProp} ) = obj.( props{jProp} );
            end

            % Get vector of Properties
            vecProperties = [];
            for jProp = 1 : length( props)
                vecProperties = [ vecProperties, obj.( props{jProp} ) ];
            end
            
            fitProperties.struct = structProperties;
            fitProperties.vec = vecProperties;

        end
        % }}}
        
        function obj = updateFittedProperties( obj, fitStruct)

            % create updated objects

        end

    end
end
