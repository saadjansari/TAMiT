% This class contains the fit machinery required to fit features to an image 
classdef FitEngine

    properties (Access = private)
        image
        feature
        parameters
    end

    methods (Access = public)
        
        % FitEngine {{{
        function obj = FitEngine( image, feature, parameters)
            % Initialize fit Engine with parameters

            obj.image = im2double(image);
            obj.feature = feature;
            obj.parameters = parameters;
        end
        % }}}

        % Set/Get methods {{{
        % Image : get and set 
        function image = GetImage(obj)
            image = obj.image;
        end
        function obj = SetImage(obj, image)
            obj.image = image;
        end

        % Feature : Get and Set
        function feature = GetFeature( obj)
            feature = obj.feature;
        end
        function obj = SetFeature(obj, feature)
            obj.feature = feature;
        end

        % Parameters : Get and Set
        function parameters = GetParameters( obj)
            parameters = obj.parameters;
        end
        function obj = SetParameters(obj, parameters)
            obj.parameters = parameters;
        end
        % }}}

        % Optimize {{{
        function obj = Optimize( obj)
            % Optimization of the feature parameters

            disp('        Optimization')

            % Local Optimization
            fitInfo = obj.OptimizeLocal;

            % Global Optimization
            fitInfo = obj.OptimizeGlobal;

            % Hyper Parameters Optimization
            fitInfo = obj.OptimizeHyperParameters( fitInfo);

        end
        % }}}

        % Methods : Optimization {{{
        
        % OptimizeLocal {{{
        function [fitInfo] = OptimizeLocal( obj)
            % Local Optimization Routine 

            if ~obj.parameters.runLocalOptimization
                fitInfo = [];
                return
            end
            
            disp('Local Optimization...')

            % Prepare for the Fit 
            [ fitProblem, fitInfo] = obj.PrepareOptimizeLocal();
            if isempty(obj.feature.featureList)
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                return
            end

            % Get total number of features to optimize 
            nFeatures = obj.feature.getSubFeatureNumber();

            for jFeature = 1 : nFeatures

                fprintf('Current Feature = %d / %d\n', jFeature, nFeatures); 

                % Solve Optimization Problem
                fitInfo{ jFeature}.fitResults = obj.SolveOptimizationProblem( fitProblem{ jFeature} );
                
                % Update feature heirarchy to ensure all dependencies are acknowledged
                obj.feature.updateSubFeatures();

            end

        end
        % }}}
        % OptimizeGlobal {{{
        function [fitInfo] = OptimizeGlobal(obj)
            % Global Optimization Routine 

            if ~obj.parameters.runGlobalOptimization
                return
            end
            
            disp('Global Optimization...')

            % Prepare for the Fit 
            [ fitProblem, fitInfo] = obj.PrepareOptimizeGlobal();
            if isempty(obj.feature.featureList)
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                return
            end

            % Solve Optimization Problem
            fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );

            % Display and Save
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end
            Cell.saveFinalFit( obj.image, obj.feature, fitInfo);

        end
        % }}}
        % OptimizeHyperParameters {{{
        function [fitInfo] = OptimizeHyperParameters( obj, fitInfoOld)
            % Optimize the hyperparameters of this fit

            % Optimize Feature Number
            fitInfo = obj.OptimizeFeatureNumber( fitInfoOld);

            % Display and Save
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end
            Cell.saveFinalFit( obj.image, obj.feature, fitInfo);

        end
        % }}}
        % OptimizeFeatureNumber {{{
        function fitInfo = OptimizeFeatureNumber( obj, fitInfoOld)
            % Optimize the number of features

            if ~obj.parameters.runFeatureNumberOptimization
                return
            end

            % check if mainFeature should be fit. If not, then return out of method
            if isempty( obj.feature.featureList)
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                fitInfo = fitInfoOld;
                return
            end

            disp('Feature Number Optimization...')

            % Book-keeping
            fitInfo.fitScope = 'globum';
            fitInfo.Org = fitInfoOld;
            fitInfo.Old = fitInfoOld;
            mainFeature.Org = obj.feature.copyDeep();
            
            % We will iteratively add and remove features to find the optimum number
            [obj, fitInfo.Old] = obj.IncreaseFeatureNumber( obj.feature, fitInfo.Old);
            [obj, fitInfo.Old] = obj.DecreaseFeatureNumber( obj.feature, fitInfo.Old);

            fprintf('- FINAL N =  %d\n', obj.feature.getSubFeatureNumber() )
            fitInfo = fitInfo.Old;
        end
        % }}}
        % IncreaseFeatureNumber {{{
        function [obj, fitInfoFinal] = IncreaseFeatureNumber( obj, feature, fitInfoOld)
            % Increase feature number until statistically insignificant
            
            continueAdd = 1;

            while continueAdd 

                % Find a residual image to use as a reference to find more features
                refImage = abs( obj.image - FitEngine.SimulateImage( fitInfoOld.fitResults.vfit, fitInfoOld) );
               
                % Create a deep copy main feature
                featureNew = feature.copyDeep();

                % Add Features
                featureNew.addSubFeatures( refImage);
                featureNew.syncFeatures();
                obj = SetFeature( obj, featureNew);
                
                fprintf('- Add N = %d ---> %d\n', feature.getSubFeatureNumber(), feature.getSubFeatureNumber()+1 )

                % Run a global fit
                [ fitProblem, fitInfo] = obj.PrepareOptimizeFeatureNumber('add');

                % Solve Optimization Problem
                fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );
                
                % Check if added feature was worth it
                im1_raw = feature.imageSim;
                im2_raw = featureNew.imageSim;
                im1 = feature.simulateFeature( size( obj.image) );
                im2 = featureNew.simulateFeature( size( obj.image) );
                numP1 = length( fitInfoOld.fitVecs.vec);
                numP2 = length( fitInfo.fitVecs.vec);
                
                % Do a statistical F test to compare the two models
                p = FitEngine.CompareModels(im1_raw, im2_raw, im1, im2, numP1, numP2, obj.image, 'f');

                % if new feature was statiscally significant, stop
                if p > obj.parameters.alpha
                    fprintf('N = %d (p = %.3f, alpha = %.3f). Keep old model.\n', feature.getSubFeatureNumber(), p, obj.parameters.alpha )
                    
                    obj = SetFeature( obj, feature);
                    fitInfoFinal = fitInfoOld;
                    continueAdd = 0;

                elseif p < obj.parameters.alpha
                    fprintf('N = %d (p = %.3f, alpha = %.3f). Keep improving model.\n',featureNew.getSubFeatureNumber(), p, obj.parameters.alpha)

                    fitInfoOld = fitInfo;
                    continueAdd = 1;
                    feature = featureNew;

                end

            end

        end
        % }}}
        % DecreaseFeatureNumber {{{
        function [obj, fitInfoFinal] = DecreaseFeatureNumber( obj, feature, fitInfoOld)
            % Increase feature number until statistically insignificant
            
            continueRemove = 1;
            while continueRemove
               
                % Find a residual image to use as a reference to find more features
                refImage = abs( obj.image - FitEngine.SimulateImage( fitInfoOld.fitResults.vfit, fitInfoOld) );
               
                % Create a deep copy main feature
                featureNew = feature.copyDeep(); 

                % Remove Features
                [~, successRemove ] = featureNew.removeSubFeatures( refImage);
                featureNew.syncFeatures();
                
                % A removal could fail if there are no features left to remove
                if successRemove

                    % Run a global fit
                    obj = SetFeature( obj, featureNew);
                    fprintf('- Rem N = %d ---> %d\n', featureNew.getSubFeatureNumber()+1, featureNew.getSubFeatureNumber() )

                    % Run a global fit
                    [ fitProblem, fitInfo] = obj.PrepareOptimizeFeatureNumber('remove');

                    % Solve Optimization Problem
                    fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );
                    
                    % Check if added feature was worth it
                    im1_raw = feature.imageSim;
                    im2_raw = featureNew.imageSim;
                    im1 = feature.simulateFeature( size( obj.image) );
                    im2 = featureNew.simulateFeature( size( obj.image) );
                    numP1 = length( fitInfoOld.fitVecs.vec);
                    numP2 = length( fitInfo.fitVecs.vec);

                    % Do a statistical F test to compare the two models
                    p = FitEngine.CompareModels(im2_raw, im1_raw, im2, im1, numP2, numP1, obj.image, 'f');

                    if p < obj.parameters.alpha
                        fprintf('N = %d (p = %.3f, alpha = %.3f). Keep old model.\n', feature.getSubFeatureNumber(), p, obj.parameters.alpha )

                        obj = SetFeature( obj, feature);
                        fitInfoFinal = fitInfoOld;
                        continueRemove = 0;

                    elseif p > obj.parameters.alpha
                        fprintf('N = %d (p = %.3f, alpha = %.3f). Keep improving model.\n',featureNew.getSubFeatureNumber(), p, obj.parameters.alpha)

                        fitInfoOld = fitInfo;
                        continueRemove = 1;
                        feature = featureNew;

                    end

                else
                    continueRemove = 0;
                    fitInfoFinal = fitInfoOld;
                    fprintf('- Rem N = %d ---> %d\n', feature.getSubFeatureNumber(), feature.getSubFeatureNumber()-1 )
                    fprintf('No more removable features\n')
                    fprintf('N = %d , Keep old model.\n', feature.getSubFeatureNumber() )
                end

            end

        end
        % }}}
        % SolveOptimizationProblem {{{
        function [fitResults, obj] = SolveOptimizationProblem(obj, problem)
            % Find the optimal solution of a problem
            
            % Run lsqnonlin fitting
            [ fitResults.vfit, ...
                fitResults.resnorm, ...
                fitResults.residual, ...
                fitResults.exitflag, ...
                fitResults.output, ~, ... 
                fitResults.jacobian] = lsqnonlin( problem); 
            fprintf('Exit Flag = %d\n', fitResults.exitflag); 
            
            % Get the errors from the jacobian
            ci = nlparci( fitResults.vfit, ...
                fitResults.residual, ...
                'jacobian', fitResults.jacobian, ...
                'alpha', 0.07); 
            fitResults.vfitError = abs( (ci(:,2)-ci(:,1) ) / 2 )';

        end
        % }}}
        
        % }}}
        
        % Methods : Prepare Optimization {{{
        
        % PrepareOptimizeLocal {{{
        function [ fitProblem, fitInformation] = PrepareOptimizeLocal( obj)
            
            % If no features are present to be fitted, assign the essentials and return out
            if isempty( obj.feature.featureList)
                fitInformation.channel = obj.parameters.channel;
                fitInformation.time = obj.parameters.time;
                fitInformation.saveDirectory = obj.parameters.saveDirectory;
                fitInformation.featureMain = obj.feature;
                fitInformation.fitScope = 'local';
                fitProblem = [];
                return
            end
        
            % Prepare local fit for all features
            nFeatures = obj.feature.getSubFeatureNumber();

            % Obtain the features to fit, their fit vector and labels 
            [ fitVec, fitLabels, fitObj] = getVecLocal( obj.feature);

            % Loop over all features
            for jFeature = 1 : nFeatures
               

                clear fitInfo
                % Get bounds of parameters to restrict the parameter space
                fitVecs = FitEngine.GetVectorBounds(obj, fitVec{jFeature},fitLabels{jFeature});


                % Scale the parameters to vary the speed of exploration in the parameter space
                if obj.parameters.fitExploreSpeed
                    speedVec = obj.getExplorationSpeedVector( fitLabels{ jFeature});
                    fitVecs.vec = fitVecs.vec ./ speedVec;
                    fitVecs.ub = fitVecs.ub ./ speedVec;
                    fitVecs.lb = fitVecs.lb ./ speedVec;
                    fitInfo.speedVec = speedVec;
                end

                % Set up parameters for fit
                fitInfo.featureMain = obj.feature;
                fitInfo.featureIndex = jFeature;
                fitInfo.fitVecs = fitVecs;
                fitInfo.mask = logical( obj.image);
                fitInfo.image = obj.image;
                fitInfo.numVoxels = size( obj.image);
                fitInfo.channel = obj.parameters.channel;
                fitInfo.time = obj.parameters.time;
                fitInfo.saveDirectory = obj.parameters.saveDirectory;
                fitInfo.alpha = 0.05;
                fitInfo.featureCurrent = fitObj{jFeature};
                fitInfo.fitScope = 'local';

                % Make the error function for optimization 
                f = obj.MakeOptimizationFcn( obj.image, fitInfo );

                % Set Optimization Options
                opts = FitEngine.SetOptimOptions( obj.parameters);
                if strcmp( obj.parameters.state, 'DEBUG')
                    plotHandle = createPlotFcnHandle( fitInfo);
                    opts = optimoptions( opts, 'OutputFcn', plotHandle);
                end

                % Create Optimization Problem
                fitProblem{ jFeature} = createOptimProblem( 'lsqnonlin', ...
                    'objective', f, ...
                    'x0', fitVecs.vec, ...
                    'ub', fitVecs.ub, ...
                    'lb', fitVecs.lb, ...
                    'options', opts);
                fitInformation{ jFeature} = fitInfo;

            end

            % Create handle for plotting function
            function stop = createPlotFcnHandle( fitInfo)
                stop = @plotFit;
                function stop = plotFit( x, optimV, state)
                    stop = Cell.plot_midFit(x, optimV, state, fitInfo);
                end
            end

        end
        % }}}
        % PrepareOptimizeGlobal {{{
        function [ fitProblem, fitInfo] = PrepareOptimizeGlobal( obj)
            
            % If no features are present to be fitted, assign the essentials and return out
            if isempty( obj.feature.featureList)
                fitInfo.channel = obj.parameters.channel;
                fitInfo.time = obj.parameters.time;
                fitInfo.saveDirectory = obj.parameters.saveDirectory;
                fitInfo.featureMain = obj.feature;
                fitInfo.fitScope = 'global';
                fitProblem = [];
                return
            end
        
            % Set Optimization Options
            opts = FitEngine.SetOptimOptions( obj.parameters);
            if strcmp( obj.parameters.state, 'DEBUG')
                opts = optimoptions( opts, 'OutputFcn', @plotFit );
            end

            % Obtain the features to fit, their fit vector and labels 
            [ fitVec, fitLabels ] = getVec( obj.feature); 
            fitObj = obj.feature; cFeature = 1; 
            fitInfo.Nnew = obj.feature.getSubFeatureNumber();

            % Get bounds of parameters to restrict the parameter space
            fitVecs = FitEngine.GetVectorBounds( obj, fitVec, fitLabels);

            % Scale the parameters to vary the speed of exploration in the parameter space
            if obj.parameters.fitExploreSpeed
                speedVec = obj.getExplorationSpeedVector( fitLabels);
                fitVecs.vec = fitVecs.vec ./ speedVec;
                fitVecs.ub = fitVecs.ub ./ speedVec;
                fitVecs.lb = fitVecs.lb ./ speedVec;
                fitInfo.speedVec = speedVec;
            end

            % Set up parameters for fit
            fitInfo.featureMain = obj.feature;
            fitInfo.featureIndex = 1;
            fitInfo.fitVecs = fitVecs;
            fitInfo.mask = logical( obj.image);
            fitInfo.image = obj.image;
            fitInfo.numVoxels = size( obj.image);
            fitInfo.channel = obj.parameters.channel;
            fitInfo.time = obj.parameters.time;
            fitInfo.saveDirectory = obj.parameters.saveDirectory;
            fitInfo.alpha = 0.05;
            fitInfo.featureCurrent = fitObj;
            fitInfo.fitScope = 'global';

            % Make the error function for optimization 
            f = obj.MakeOptimizationFcn( obj.image, fitInfo );

            % Create Optimization Problem
            fitProblem = createOptimProblem( 'lsqnonlin', ...
                'objective', f, ...
                'x0', fitVecs.vec, ...
                'ub', fitVecs.ub, ...
                'lb', fitVecs.lb, ...
                'options', opts);

            % Create handle for plotting function
            function stop = plotFit( x, optimV, state)
                stop = Cell.plot_midFit(x, optimV, state, fitInfo);
            end

        end
        % }}}
        % PrepareOptimizeFeatureNumber {{{
        function [ fitProblem, fitInfo] = PrepareOptimizeFeatureNumber( obj, routine)
            % routine = 'add' or 'remove'
            
            [ fitProblem, fitInfo] = obj.PrepareOptimizeGlobal;

            switch routine
                case 'add'
                    fitInfo.fitScope = 'globum_add';
                    fitInfo.Nold = fitInfo.Nnew - 1;

                case 'remove'
                    fitInfo.fitScope = 'globum_remove';
                    fitInfo.Nold = fitInfo.Nnew + 1;
            end

            opts = fitProblem.options;
            if strcmp( obj.parameters.state, 'DEBUG')
                opts = optimoptions( opts, 'OutputFcn', @plotFit );
            end

            % Make the error function for optimization 
            f = obj.MakeOptimizationFcn( obj.image, fitInfo );

            % Create Optimization Problem
            fitProblem = createOptimProblem( 'lsqnonlin', ...
                'objective', f, ...
                'x0', fitInfo.fitVecs.vec, ...
                'ub', fitInfo.fitVecs.ub, ...
                'lb', fitInfo.fitVecs.lb, ...
                'options', opts);

            % Create handle for plotting function
            function stop = plotFit( x, optimV, state)
                stop = Cell.plot_midFit(x, optimV, state, fitInfo);
            end

        end
        % }}}
        
        % }}}
        
    end

    methods (Static = true)

        % SetOptimOptions {{{
        function opts = SetOptimOptions( config)
            % config.state      = DEBUG 
            %                   = RELEASE
            % config.parallel   = true/false

            % Set default options
            opts = optimoptions( @lsqnonlin, ...
                                'MaxFunEvals', 2000, ...
                                'OptimalityTolerance', 1e-12, ...
                                'MaxIter', 20, ...
                                'TolFun', 1e-9, ...
                                'FiniteDifferenceStepSize', 1e-2, ...
                                'FiniteDifferenceType', 'central', ...
                                'StepTolerance', 1e-5);

            % Set configurable options
            switch config.state
                case 'RELEASE'
                    opts = optimoptions( opts, ...
                                        'display', 'off' );

                case 'DEBUG'
                    opts = optimoptions( opts, ...
                                        'display', 'iter' );
            end
            
            % Set Parallel Optimization
            if config.useParallel
                opts = optimoptions( opts, 'UseParallel', true); 
            end
            
        end
        % }}}

        % SimulateImage{{{
        function imageOut = SimulateImage( vec, fitInfo)
            % Simulates an image with the feature from fitInfo
            
            % Apply the speedExploration vector to transform to real parameters
            if isfield( fitInfo, 'speedVec')
                vec = vec .* fitInfo.speedVec;
            end

            % Get feature to simulate
            feature = fitInfo.featureCurrent;
            featureID = feature.ID;

            % Update feature with vector parameters
            feature.absorbVec( vec, fitInfo.fitVecs.labels );
            
            % Simulate the feature
            imageOut = fitInfo.featureMain.simulateAll( 0*fitInfo.image, featureID);

        end
        % }}}
        
        % CompareModels {{{
        function vargout = CompareModels(im1_raw, im2_raw, im1, im2, numP1, numP2, imOrig, test)
            
            % Assumes that feature images are without the background
            
            if ~strcmp( test, 'f') && ~strcmp( test, 'aic') && ~strcmp( test, 'bic')
                error('compareModels: unknown test (must be ''f'', ''aic'' or ''bic''.')
            end
            
            numObs1 = sum( im1(:) > 0);
            numObs2 = sum( im2(:) > 0);
            
            dof1 = numObs1-numP1;
            dof2 = numObs2-numP2;
            
            resid1 = (imOrig - im1).^2;
            resid2 = (imOrig - im2).^2;
            res1 = sum( resid1(:) );
            res2 = sum( resid2(:) );
            
            % residual at feature regions
            % find feature region
            im1_raw( im1_raw < 0.001*max(im1_raw(:) ) ) = 0;
            im2_raw( im2_raw < 0.001*max(im2_raw(:) ) ) = 0;
            ROI = logical( im1_raw > 0 | im2_raw > 0 );
            resR1 = sum( resid1( ROI) );
            resR2 = sum( resid2( ROI) );

            switch test
                case 'f'
%                     p1 = fcdf( res2/res1 , dof2, dof1);
                    
                    dof1 = sum( ROI(:) )-numP1;
                    dof2 = sum( ROI(:) )-numP2;
                    p2 = fcdf( resR2/resR1 , dof2, dof1);
                    
                    vargout = p2;
                    
                case 'aic'
                    
                case 'bic'
                    
            end
            
            
        end
        % }}}

        % MakeOptimizationFcn {{{
        function errVal = MakeOptimizationFcn( ImageFit, fitInfo)

            errVal = @OptimizationFcn;

            function err = OptimizationFcn( p)
                err = FitEngine.SimulateImage( p, fitInfo) - ImageFit;
                err = err(:);
            end

        end
        % }}}
        
        % GetVectorBounds {{{
        function fitVecs = GetVectorBounds( obj, vec, vecLabels)
            % Get upper and lower bounds for the parameters in the fit vector
            
            image = obj.GetImage();
            ub = vec; lb = vec;

            % Spatial 
            minVox = 1;
            maxVox = size( image);

            % Background
            estBkg = median( image( image > 0) );
            maxBkg = max( image(:) );
            minBkg = min( image(:) );

            % Amplitude
            maxAmp = max( image(:) );
            minAmp = 0; 

            % Sigma
            minSig = [1.2, 1.2, 1.0];
            maxSig = [2.0, 2.0, 2.0];

            % Find the correct label in vecLabels, and place the correct bounds in the correct places

            % Find index of parameters 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'amplitude') ) );
            idxSig = find( ~cellfun( @isempty, strfind( vecLabels, 'sigma') ) );
            idxP0 = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            idxP1 = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            idxP = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            idxBkg = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );

            % Store upper and lower bounds correctly
            if ~isempty( idxAmp), 
                ub( idxAmp) = maxAmp;
                lb( idxAmp) = minAmp; end
            if ~isempty( idxSig), 
                nF = length( idxSig)/3;
                ub( idxSig) = repmat( maxSig, 1, nF);
                lb( idxSig) = repmat( minSig, 1, nF); end
            if ~isempty( idxP0), 
                nF = length( idxP0)/3;
                ub( idxP0 ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP0 ) = minVox; end
            if ~isempty( idxP1), 
                nF = length( idxP1)/3;
                ub( idxP1 ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP1 ) = minVox; end
            if ~isempty( idxP), 
                nF = length( idxP)/3;
                ub( idxP ) = repmat( [ maxVox(2) maxVox(1) maxVox(3) ], 1, nF);
                lb( idxP ) = minVox; end
            if ~isempty( idxBkg), 
                ub( idxBkg ) = maxBkg;
                lb( idxBkg ) = minBkg; end

            if any( lb > ub) || any(ub < lb) || any( vec < lb) || any(vec > ub)
                badIdx = unique([ find( lb > ub),  find(ub < lb) , find( vec < lb) , find(vec > ub)]);
                disp('Bad prop labels')
                badProps = vecLabels( badIdx)
                disp(vecLabels)
                disp( ub )
                disp( vec)
                disp( lb )
                error('getUpperLowerBounds : bounds are wrong')
            end

            fitVecs.vec = vec;
            fitVecs.ub = ub;
            fitVecs.lb = lb;
            fitVecs.labels = vecLabels;

        end
        % }}}
        
        % getExplorationSpeedVector {{{
        function speedVec = getExplorationSpeedVector( vecLabels)
            % Creates a weighing vector to allow a user to assign different importance to different kinds of parameters. This will infact allow the fitting optimization engine to explore the parameters space at different speeds
            
            % Exlporation Speed : unassigned speeds are kept at 1.0
%             speedAmp = 100;
%             speedBkg = 10;
%             speedSigma = 1;
%             speedPos = 10;
            speedAmp = 1000;
            speedBkg = 10;
            speedSigma = 1;
            speedPos = 1;
            speedVec = ones( size(vecLabels) );

            % Find the index of these speeds
            % Amplitude
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'amplitude') ) );
            speedVec( idxAmp) = speedAmp;

            % Background 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );
            speedVec( idxAmp) = speedBkg;

            % Sigma 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'sigma') ) );
            speedVec( idxAmp) = speedSigma;
            
            % Position
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            speedVec( idxAmp) = speedPos;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            speedVec( idxAmp) = speedPos;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            speedVec( idxAmp) = speedPos;

        end
        % }}}
        
    end

end
