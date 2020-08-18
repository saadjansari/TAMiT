classdef MonopolarAster < OrganizerMaster
% A MonopolarAster is a higher level feature that sits inside a monopolar cell. It is composed of 1 AsterMT 
    properties
    end

    methods (Access = public)
        
        % MonopolarAster {{{
        function obj = MonopolarAster( dim, image, aster, props2Fit)

            % ensure that there is only 1 element in featureList and that is of type AsterMT
            if length( aster) ~= 1 || ~strcmp( aster{1}.type, 'AsterMT')
                error('MonopolarAster: aster must be of type ''AsterMT'' ') 
            end

            obj = obj@OrganizerMaster( dim, image, aster, props2Fit, 'MonopolarAster');

        end
        % }}}
        
        % getVec {{{
        function [vec, vecLabels,ub,lb] = getVec( obj, props)

            % Get general environmental properties
            [vecE, vecLabelsE] = obj.getVecEnvironment();
            ubE = vecE; lbE = vecE;

            % get vectors from features
            props.spb = {'position', 'amplitude', 'sigma'};
            props.mt = {'endPosition','length', 'theta','amplitude', 'sigma'};

            % get vectors from aster. 
            [vec, vecLabels, ub, lb] = getVec( obj.featureList{1}, props.spb, props.mt);
%             vec = [vecE, vec]; ub = [ubE, ub]; lb = [lbE, lb]; 
%             vecLabels = { vecLabelsE{:}, vecLabels{:} };
        end
        % }}}

        % absorbVec {{{
        function obj = absorbVec( obj, vec, vecLabels)
            
            % Absorb Environmental parameters
            obj.absorbVecEnvironment( vec, vecLabels);

            % Aster Vector
            obj.featureList{ 1}.absorbVec( vec, vecLabels );
            
        end
        % }}}

        % getVecLocal {{{
        function [ vecList, vecLabelsList, objList, ubList, lbList] = getVecLocal( obj)

            % I can have a list of vectors and a list of vectorLabels, and a list of objects
            % Order of features : 
            %   Spindle Microtubule (complete optimization)
            %   Astral Microtubules (end optimization)

            vecList = {};
            vecLabelsList = {};
            objList = {};
            ubList = {}; lbList = {};

            % Define properties to get
            props.spb = {'position', 'amplitude', 'sigma'};
%             props.mt = {'endPosition', 'amplitude', 'sigma'};
            props.mt = {'length','theta', 'amplitude', 'sigma'};
            % SPB vector
            try
                [ vecList{1}, vecLabelsList{1}, ubList{1}, lbList{1}] = getVec( obj.featureList{1}.featureList{1}, props.spb);
            catch
                [ vecList{1}, vecLabelsList{1}] = getVec( obj.featureList{1}.featureList{1}, props.spb);
            end
            objList{1} = obj.featureList{1}.featureList{1};
            numFeatures = 1;

            % Microtubule vectors : Aster
            for jmt = 1 : obj.featureList{1}.numFeatures-1
                lab = numFeatures+jmt;
                try
                    [ vecList{ lab} , vecLabelsList{ lab}, ubList{lab}, lbList{lab}] = getVec( obj.featureList{1}.featureList{1+jmt}, props.mt );
                catch
                    [ vecList{ lab} , vecLabelsList{ lab}] = getVec( obj.featureList{1}.featureList{1+jmt}, props.mt );
                end
                objList{ lab} = obj.featureList{1}.featureList{1+jmt};
            end
            numFeatures = numFeatures + obj.featureList{1}.numFeatures-1;

            % Prepend environmental parameters to each vec and vecLabel
%             [vecE, vecLabelsE] = obj.getVecEnvironment();
%             
%             for jFeat = 1 : numFeatures
%                 vecList{jFeat} = [ vecE , [vecList{jFeat}(:)]' ];
%                 vecLabelsList{jFeat} = { vecLabelsE{:} , vecLabelsList{jFeat}{:} };
%             end

        end
        % }}}

        % addSubFeatures {{{
        function [obj, successAdd] = addSubFeatures( obj, Image2Find)
            % Searches for missing features and adds them to the featureList
            
            % For a monopolar aster, we can only add more astral microtubules : we will find a possible microtubule to add 

            % How to find a possible microtubule:
            %   1. Angular Sweep (not clear how to ensure a line is found if there are no peaks)
            %   2. Pick the location of the highest residual and draw a line from both pole connecting to the max residual point. Find the mean intensity of each line. Subtract the mean from the actual voxel values along the line to find a total residual. Pick the pole with the minimum residual.

            props = Cell.GetFeatureProps();
            props2Fit = {'endPosition', 'theta','amplitude', 'sigma'};
            display = props.line.graphics.green;
            if obj.dim==3, sigma=[1.2 1.2 1.0]; elseif obj.dim==2, sigma=[1.2 1.2]; end

            % ask subfeatures to find missing features
            Image2Find( Image2Find < 0) = 0;
            [featureKeep, succ] = obj.featureList{1}.findBestMissingFeature( obj.image, Image2Find, [], []);
            
            if ~succ
                successAdd = 0; return
            end
            
            % Create feature object 
            amp = featureKeep.amplitude;
            feature = Line( featureKeep.startPosition, featureKeep.endPosition, amp, sigma, obj.dim, props2Fit, display);
            feature.repr = 'spherical';
            if feature.length < 7
                successAdd = 0;
                return
            end
            
            % Add feature to the correct subfeature 
            idxAdd = obj.featureList{1}.addFeatureToList( feature);
            successAdd = 1;
            
        end
        % }}}

        % removeSubFeatures {{{
        function [ obj, successRemove] = removeSubFeatures( obj, Image2Find)
            % Searches for redundant features and removes them from the featureList
            
            % Ask aster for least prominent feature
            [feat, successRemove] = obj.featureList{1}.findLeastPromFeature();
            
            if successRemove
                % Remove the worst microtubule
                obj.featureList{1}.removeFeatureFromList( feat.idx);
            end
            
        end
        % }}}

        % updateSubFeatures {{{
        function obj = updateSubFeatures(obj)
            
            % ensure aster pole is same as start positions of all microtubules
            aster = obj.featureList{1};
            pole = aster.featureList{1};
            pos = pole.position;

            for jmt = aster.numFeatures-1
                aster.featureList{jmt+1}.startPosition = pos;
            end

        end
        % }}}
        
        function obj = finalizeAddedFeatures(obj)
            
            % Display properties for lines in asters
            display = {'Color', [1 0 1], 'LineWidth', 3};
            
            for jL = 2 : length(obj.featureList{1} )
                obj.featureList{1}.featureList{jL}.display = display;
            end
                
        end
        
    end

    methods( Static = true )

        % loadFromStruct {{{
        function obj = loadFromStruct( S)
            
            if ~isfield( S, 'type') || ~strcmp( S.type, 'MonopolarAster')
                error('incorrect type')
            end

            % Load all the features from their structures recursively
            featureList{ 1} = AsterMT.loadFromStruct( S.featureList{ 1} ); 

            obj = MonopolarAster( S.dim, S.image, featureList, S.props2Fit);
%             obj = obj@Organizer( S.dim, featureList, S.type);
            obj.findEnvironmentalConditions();
            %obj.syncFeaturesWithMap();

        end
        % }}}

    end

end
