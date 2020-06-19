classdef InterphaseCell < Cell
% This is a specialized cell of superclass Cell that is used for an interphase cell.
    properties
        phase = 'interphase'
    end

    methods
        
        % InterphaseCell {{{
        function obj = InterphaseCell( imageData, featuresInChannels, channelsToFit, params, varargin)
        % InterphaseCell: constructor function for InterphaseCell object     

            obj = obj@Cell( imageData, featuresInChannels, channelsToFit, params, varargin{:}, 'Type', 'Interphase');

        end
        % }}}

        % EstimateFeatures {{{
        function obj = EstimateFeatures( obj, estimationImage, cTime, cChannel, idxChannel, timeReverse)
        % findFeatures : estimates and finds the features 
            
            % Get feature type
            currentFeature  = obj.featuresInChannels{ idxChannel};

            % Get Start time 
            lifetime = obj.imageData.GetLifetime;
            startTime = lifetime(1);

            % Novel Estimation for first frame
            if cTime == startTime
                obj.featureList{ idxChannel, cTime} = obj.EstimateFeaturesNovel( currentFeature, estimationImage);
            % Propagate old feature for later frames
            else 
                obj = obj.PropagateOldFeature( idxChannel, cChannel, cTime);
            end

        end
        % }}}

        % EstimateFeaturesNovel {{{
        function feature = EstimateFeaturesNovel( obj, currentFeature, image)
            disp('- DeNovo') 

            switch currentFeature
                case 'Microtubule'
                    feature = IMTBank.findIntBank( im2double(image), obj.params.estimate.interphase, obj.featureProps.intBank);
                    %feature = InterphaseCell.findFeaturesDeNovo_MT( image, obj.params.estimate.interphase);

                %case 'Kinetochore'
                    %feature = MonopolarCell.findFeaturesDeNovo_KC( image, obj.params.estimate.kcbank);

                %case 'Cut7'
                    %feature = MonopolarCell.findFeaturesDeNovo_Cut7( image, obj.params.estimate.cut7dist);
            end

        end
        % }}}

    end

    methods ( Static = true, Access = public )

        % Microtubules {{{
        
        % findFeaturesDeNovo_MT {{{
        function iMTBankObj = findFeaturesDeNovo_MT( imageIn, params)
            % Interphase Cell:
            % Find special interphase asters (MTOC + 2 curves in an opposite direction)

            dim = length( size(imageIn) );
            imageIn = im2double( imageIn);
            bkg = median( imageIn( imageIn(:) > 0) );
            
            if dim == 3
                sigma=[1.2 1.2 1.0];
                props.iMT = {'cX', 'cY', 'cZ', 'amplitude', 'sigma'};
            elseif dim == 2
                sigma=[1.2 1.2];
                props.iMT = {'cX', 'cY', 'amplitude', 'sigma'};
            end
            props.MTOC= {'position', 'amplitude', 'sigma'};
            props.iMTBank = {'background'};
            
            display.MTOC= {'Color', [0.7 0 0.7] , 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2};
            display.iMT= {'Color', [1 0 1] , 'LineWidth', 2};
            
            % Search for asters
            iAsters = IMTBank.findAstersInterphase( imageIn, params, sigma, props, display);

            % Store AsterMT stored in a iMTBank 
            iMTBankObj = IMTBank( dim, imageIn, iAsters, props.iMTBank);
            iMTBankObj.findEnvironmentalConditions();
            
        end
        % }}}

        % }}}

    end

end
