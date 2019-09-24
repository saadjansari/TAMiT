classdef Cut7Distribution < Organizer 
% This is the main organizer of the cut7 distribution. 
    properties
        image
        spindlePositionStart
        spindlePositionEnd
    end

    methods (Access = public)

        % Cut7Distribution {{{
        function obj = Cut7Distribution( image, featureList)

            if nargin < 1 
                error( 'Cut7Distribution: number of inputs is inconsistent')
            elseif nargin < 1
                warning('Cut7Distirbution: initialized as an empty container')
            end

            dim = length( size( image) );
            obj = obj@Organizer( dim, featureList, 'Cut7Distribution');
            obj.image = image;

            % initialize featureMap and assign IDs
            %obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

        end
        % }}}
        
        % saveAsStruct {{{
        function S = saveAsStruct( obj)

            S.type = obj.type;
            S.image = uint16( obj.image);

            for jFeat = 1 : length( obj.featureList)
                S.featureList{ jFeat} = saveAsStruct( obj.featureList{jFeat} );
            end

            S.spindlePositionStart = obj.spindlePositionStart;
            S.spindlePositionEnd = obj.spindlePositionEnd;

        end
        % }}}
        
        % harvestSpindleInfo {{{
        function harvestSpindleInfo( obj, spindleObj)
            % gets spindle information from spindleObj and stores it
            
            obj.spindlePositionStart = spindleObj.featureList{1}.startPosition;
            obj.spindlePositionEnd = spindleObj.featureList{1}.endPosition;
            disp( 'Harvested spindle information')

        end
        % }}}

    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'Cut7Distribution')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            %for jFeat = 1 : length( S.featureList)
                %featureList{ jFeat} = Blob.loadFromStruct( S.featureList{ jFeat} ); 
            %end

            obj = Cut7Distribution( S.image, {});
            obj.spindlePositionStart = S.spindlePositionStart;
            obj.spindlePositionEnd = S.spindlePositionEnd;

        end
        % }}}

    end

end
