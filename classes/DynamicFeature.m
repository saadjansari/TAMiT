classdef DynamicFeature
    % This defines a dynamic feature that exists over multiple time frames

    % properties {{{
    properties
        featureType
        featureProps
        startTime
        endTime
        sourceData
        displayProps
    end
    % }}}

    % Main Methods {{{
    methods
    
        % Initialization
        % Dynamic Feature {{{
        function obj = DynamicFeature( featureType, startTime, endTime, sourceData, featureProps, displayProps)
       
            % FeatureType specific assignment {{{
            % acceptable feature types: spots, lines (straight or curved)
            obj.featureType = featureType;
            if strcmp( featureType, 'spots')
                
                % check if featureProps has the relevant information
                fields = {'coord', 'amp', 'std'};
                for jfield = 1 : length( fields)
                    if ~isfield( featureProps, fields{jfield} )
                        error( ['DynamicFeature: the ''spots'' feature must have the field : ' fields{jfield} ] )
                    end
                end

                obj.featureProps.coord = featureProps.coord;
                obj.featureProps.amp= featureProps.amp;
                obj.featureProps.std= featureProps.std;

            elseif strcmp( featureType, 'lines')

                % check if featureProps has the relevant information
                fields = {'polyCoef', 'amp', 'std'};
                for jfield = 1 : length( fields)
                    if ~isfield( featureProps, fields{jfield} )
                        error( ['DynamicFeature: the ''lines'' feature must have the field : ' fields{jfield} ] )
                    end
                end

                obj.featureProps.polyCoef= featureProps.polyCoef;
                obj.featureProps.amp= featureProps.amp;
                obj.featureProps.std= featureProps.std;

            end
            % }}}

            % Times
            obj.startTime = startTime;
            obj.endTime = endTime;

            % Source Data
            obj.sourceData = sourceData;

            % Display Props
            obj.displayProps = displayProps;

        end
        % }}}

    end
    % }}}

    % Spot Methods {{{
    methods


    end
    % }}}

    % Line Methods {{{
    methods

        % Line length {{{
        function obj = findLineLength( obj)


        end
        % }}}

        % Line Orientation {{{
        function obj = findLineOrientation( obj)

        end
        % }}}

    end
    % }}}

end
