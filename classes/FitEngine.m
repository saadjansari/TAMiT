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

        % Methods : Optimization {{{
        
        % Optimize {{{
        function obj = Optimize( obj)
            % Optimization of the feature parameters

            disp('        Optimization')

            fitInfo = [];
            if obj.parameters.fit2DFirst == 1
                [obj,fitInfo] = obj.OptimizeProjection2D();
            end

            % Environment Optimization
            fitInfo = obj.OptimizeEnvironment(fitInfo); fitInfo = [];
            
            % Local Optimization
            fitInfo = obj.OptimizeLocal(fitInfo);
            
            % Global Optimization
            fitInfo = obj.OptimizeGlobal(fitInfo);
            
            % Hyper Parameters Optimization
            fitInfo = obj.OptimizeHyperParameters( fitInfo);
            
            % Save final optimized data
            fitInfo.fitScope = 'final';
            obj.feature = fitInfo.featureMain;
            obj.feature.finalizeAddedFeatures();
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end
            Cell.saveFinalFit( obj.image, obj.feature, fitInfo);


        end
        % }}}
        % OptimizeLocal {{{
        function [fitInfo] = OptimizeLocal( obj, fitInfo)
            % Local Optimization Routine 

            if ~obj.parameters.runLocalOptimization
                return
            end
            if isempty(obj.feature.featureList)
                fitInfo.channel = obj.parameters.channel;
                fitInfo.time = obj.parameters.time;
                fitInfo.saveDirectory = obj.parameters.saveDirectory;
                fitInfo.featureMain = obj.feature;
                fitInfo.fitScope = 'local';
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                return
            end
            
            disp('Local Optimization...')


            % Get total number of features to optimize 
            nFeatures = obj.feature.getSubFeatureNumber();

            for jFeature = 1 : nFeatures

                fprintf('Current Feature = %d / %d\n', jFeature, nFeatures);

                % Prepare for optimization
                obj.feature.forceInsideMask();
                [ fitProblem{jFeature}, fitInfo{jFeature}] = obj.PrepareOptimizeLocal( jFeature);

                % Solve Optimization Problem
%                 if jFeature == 4
%                     stoph = 1;
%                 end
                fitInfo{ jFeature}.fitResults = obj.SolveOptimizationProblem( fitProblem{ jFeature} );
                obj.updateFinalFeature( fitInfo{jFeature} );
                
                % Update feature heirarchy to ensure all dependencies are acknowledged
                obj.feature.updateSubFeatures();

            end

        end
        % }}}
        % OptimizeGlobal {{{
        function [fitInfo] = OptimizeGlobal(obj, fitInfo)
            % Global Optimization Routine 

            if ~obj.parameters.runGlobalOptimization
                return
            end
            
            disp('Global Optimization...')

            % Prepare for the Fit
            obj.feature.forceInsideMask();
            [ fitProblem, fitInfo] = obj.PrepareOptimizeGlobal();
            if isempty(obj.feature.featureList)
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                return
            end

            % Solve Optimization Problem
            fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );
            obj.updateFinalFeature( fitInfo);
            fitInfo.fitInfoOld = fitInfo;
            fitInfo.fitVecOld = fitInfo.fitVecs.vec;

            % Display and Save
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end
            %Cell.saveFinalFit( obj.image, obj.feature, fitInfo);

        end
        % }}}
        % OptimizeEnvironment {{{
        function [fitInfo] = OptimizeEnvironment(obj, fitInfo)
            % Global Optimization Routine 

            if ~obj.parameters.runGlobalOptimization
                return
            end
            
            disp('Environment Optimization...')

            % Prepare for the Fit
            obj.feature.forceInsideMask();
            [ fitProblem, fitInfo] = obj.PrepareOptimizeEnvironment();
            if isempty(obj.feature.featureList)
                fprintf('Skipping fitting for %s\n', obj.feature.type);  
                return
            end

            % Solve Optimization Problem
            fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );
            obj.feature.absorbVecEnvironment( fitInfo.fitResults.vfit.*fitInfo.speedVec, fitInfo.fitVecs.labels);
            fitInfo.fitInfoOld = fitInfo;
            fitInfo.fitVecOld = fitInfo.fitVecs.vec;

            % Display and Save
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end

        end
        % }}}
        % OptimizeHyperParameters {{{
        function [fitInfo] = OptimizeHyperParameters( obj, fitInfoOld)
            % Optimize the hyperparameters of this fit

            % Optimize Feature Number
            [obj,fitInfo] = obj.OptimizeFeatureNumber( fitInfoOld);
            fitInfo.fitInfoOld = fitInfoOld;
            fitInfo.fitVecOld = fitInfoOld.fitResults.vfit;

            % Display and Save
%             if obj.parameters.display
%                 Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
%             end
            %Cell.saveFinalFit( obj.image, obj.feature, fitInfo);

        end
        % }}}
        % OptimizeFeatureNumber {{{
        function [obj,fitInfo] = OptimizeFeatureNumber( obj, fitInfoOld)
            % Optimize the number of features

            if ~obj.parameters.runFeatureNumberOptimization
                fitInfo = fitInfoOld;
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
%             mainFeature.Org = obj.feature.copyDeep();
            
            % We will iteratively add and remove basic features to find the optimum number
            disp('Feature Number Optimization...')
            [obj, fitInfo.Old] = obj.IncreaseFeatureNumber( obj.feature, fitInfo.Old, 'Basic');
            [obj, fitInfo.Old] = obj.DecreaseFeatureNumber( obj.feature, fitInfo.Old, 'Basic');
            
            fprintf('- FINAL N =  %d\n', obj.feature.getSubFeatureNumber() )                                   
            fitInfo = fitInfo.Old;
        end
        % }}}
        % IncreaseFeatureNumber {{{
        function [obj, fitInfoFinal] = IncreaseFeatureNumber( obj, feature, fitInfoOld, featureType)
            % Increase feature number until statistically insignificant
            
            continueAdd = 1;

            while continueAdd 

                % Find a residual image to use as a reference to find more features
                refImage = abs( obj.image - FitEngine.SimulateImage( fitInfoOld.fitResults.vfit, fitInfoOld) );
               
                % Create a deep copy main feature
                featureNew = feature.copyDeep();

                % Add Features
                switch featureType
                    case 'Basic'
                        [~, successAdd ] = featureNew.addSubFeatures( refImage);
                    case 'Organizer'
                        [~, successAdd ] = featureNew.addOrganizers( refImage);
                end
                
                featureNew.syncFeatures();
                
                % An addition could fail if there are no features that
                % could be added
                if successAdd
                    
                    featureNew.forceInsideMask();
                    obj = SetFeature( obj, featureNew);

                    fprintf('- Add N = %d ---> %d\n', feature.getSubFeatureNumber(), featureNew.getSubFeatureNumber() )

                    % Run a global fit
                    [ fitProblem, fitInfo] = obj.PrepareOptimizeFeatureNumber('add');

                    % Solve Optimization Problem
                    fitInfo.fitResults = obj.SolveOptimizationProblem( fitProblem );

                    % Check if added feature was worth it
                    im1_raw = feature.simulateFeature( size(obj.image));
                    im2_raw = featureNew.simulateFeature( size(obj.image));
                    im1 = feature.simulateAll( 0*obj.image, feature.ID );
                    im2 = featureNew.simulateAll( 0*obj.image, featureNew.ID );
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
                    
                else
                    continueAdd = 0;
                    fitInfoFinal = fitInfoOld;
                    fprintf('- Add N = %d +\n', feature.getSubFeatureNumber() )
                    fprintf('No more features to add\n')
                    fprintf('N = %d , Keep old model.\n', feature.getSubFeatureNumber() )
                end

            end

        end
        % }}}
        % DecreaseFeatureNumber {{{
        function [obj, fitInfoFinal] = DecreaseFeatureNumber( obj, feature, fitInfoOld, featureType)
            % Increase feature number until statistically insignificant
            
            continueRemove = 1;
            while continueRemove
               
                % Find a residual image 
                refImage = abs( obj.image - FitEngine.SimulateImage( fitInfoOld.fitResults.vfit, fitInfoOld) );
               
                % Create a deep copy main feature
                featureNew = feature.copyDeep(); 

                % Remove Features
                switch featureType
                    case 'Basic'
                        [~, successRemove ] = featureNew.removeSubFeatures( refImage);
                    case 'Organizer'
                        [~, successRemove ] = featureNew.removeOrganizers( refImage);
                end
                featureNew.updateFeatureMap();
                
                if featureNew.getSubFeatureNumber() == 0
                    successRemove = 0;
                end
                
                % A removal could fail if there are no features left to remove
                if successRemove

                    % Run a global fit
                    featureNew.forceInsideMask();
                    obj = SetFeature( obj, featureNew);
                    fprintf('- Rem N = %d ---> %d\n', feature.getSubFeatureNumber(), featureNew.getSubFeatureNumber() )

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
                    fprintf('- Rem N = %d -\n', feature.getSubFeatureNumber() )
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
        % OptimizeProjection2D {{{
        function [obj,fitInfo] = OptimizeProjection2D( obj)
           % fit 2D projection of the features
           
            % Only proceed if dimensionality is 3
            dim = length( size(obj.image ));
            if dim == 2
                disp('Skipping 2D projection fitting: current dimensionality is 2')
                return
            end
            
            % Get 2D image
            Image2D = max( obj.image, [], 3);
            
            % Get 2D projection of features
            feature2D = obj.feature.copyDeep();
            feature2D = feature2D.GetProjection2D();
            feature2D.fillParams( size(Image2D) );
            
            % Initialize 2D fit engine
            params = obj.parameters;
            params.runLocalOptimization = 0;
            params.runGlobalOptimization = 1;
            params.runFeatureNumberOptimization = 1;
            fitEngine2 = FitEngine( Image2D, feature2D, params);
            
            % Fit global 2D projection of feature
            fitInfo = fitEngine2.OptimizeGlobal();
            obj.parameters.runLocalOptimization = 0;
            
%             if params.runFeatureNumberOptimization
%                 fitInfo = fitEngine2.OptimizeHyperParameters(fitInfo);
% %                 obj.feature = fitEngine2.feature;
%                 obj.parameters.runFeatureNumberOptimization = 0;
%             end
            
            % Update 3D features using fitted 2D information
            obj.feature.Update3DFrom2D(fitInfo.featureCurrent);
            obj.feature.syncFeatures();
            
        end
        % }}}
        
        % }}}
        
        % Methods : Prepare Optimization {{{
        
        % PrepareOptimizeLocal {{{
        function [ fitProblem, fitInfo] = PrepareOptimizeLocal( obj, jFeature)
            
            % Prepare local fit for all features
            nFeatures = obj.feature.getSubFeatureNumber();

            % Obtain the features to fit, their fit vector and labels 
            try
                [ fitVec, fitLabels, fitObj, ubList, lbList] = getVecLocal( obj.feature);
            catch
                [ fitVec, fitLabels, fitObj] = getVecLocal( obj.feature);
                ubList = cell(1, length(fitVec)); lbList = cell(1, length(fitVec));
            end
            % Get bounds of parameters to restrict the parameter space
            fitVecs = FitEngine.GetVectorBounds(obj, fitVec{jFeature},fitLabels{jFeature}, ubList{jFeature}, lbList{jFeature});

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
            fitInfo.featureCurrent.fillParams( size( fitInfo.image));
            fitInfo.timeReversal = obj.parameters.timeReversal;

            % Make the error function for optimization 
            f = obj.MakeOptimizationFcn( obj.image, fitInfo );

            % Set Optimization Options
            opts = FitEngine.SetOptimOptions( obj.parameters);
            if strcmp( obj.parameters.state, 'DEBUG')
                plotHandle = createPlotFcnHandle( fitInfo);
                fillParamHandle = createFillParamsFcnHandle(fitInfo);
                opts = optimoptions( opts, 'OutputFcn', {plotHandle,fillParamHandle});
            else
                fillParamHandle = createFillParamsFcnHandle(fitInfo);
                opts = optimoptions( opts, 'OutputFcn', fillParamHandle);
            end

            % Create Optimization Problem
            fitProblem = createOptimProblem( 'lsqnonlin', ...
                'objective', f, ...
                'x0', fitVecs.vec, ...
                'ub', fitVecs.ub, ...
                'lb', fitVecs.lb, ...
                'options', opts);

            % Create handle for plotting function
            function stop = createPlotFcnHandle( fitInfo)
                stop = @plotFit;
                function stop = plotFit( x, optimV, state)
                    stop = Cell.plot_midFit(x, optimV, state, fitInfo);
                end
            end
            function stop = createFillParamsFcnHandle( fitInfo)
                stop = @fillParams;
                function stop = fillParams( x, optimV, state)
                    stop = false;
                    fitInfo.featureCurrent.fillParams( size( fitInfo.image));
                end
            end

        end
        % }}}
        % PrepareOptimizeEnvironment {{{
        function [ fitProblem, fitInfo] = PrepareOptimizeEnvironment( obj)
            
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
            %if strcmp( obj.parameters.state, 'DEBUG')
                %opts = optimoptions( opts, 'OutputFcn', {@plotFit, @fillParams} );
            %else
                %opts = optimoptions( opts, 'OutputFcn', @fillParams );
            %end

            % Obtain the features to fit, their fit vector and labels 
            % Get bounds of parameters to restrict the parameter space
            [fitVec, fitLabels] = obj.feature.getVecEnvironment();
            ub = max( obj.image(:)); lb = 0;
            fitVecs.vec = fitVec; fitVecs.labels = fitLabels; fitVecs.ub = ub; fitVecs.lb = lb;

            fitObj = obj.feature; cFeature = 1; 
            fitInfo.Nnew = obj.feature.getSubFeatureNumber(); 
            
            % Scale the parameters to vary the speed of exploration in the parameter space
            if obj.parameters.fitExploreSpeed
                speedVec = obj.getExplorationSpeedVector( fitLabels);
                fitVecs.vec = fitVecs.vec ./ speedVec;
                fitVecs.ub = fitVecs.ub ./ speedVec;
                fitVecs.lb = fitVecs.lb ./ speedVec;
                fitInfo.speedVec = speedVec;
            end

            % Set up parameters for fit
            fitInfo.optimizeEnvOnly = 1;
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
            fitInfo.featureCurrent.fillParams( size( fitInfo.image));
            fitInfo.timeReversal = obj.parameters.timeReversal;

            % Make the error function for optimization 
            f = obj.MakeOptimizationFcn( obj.image, fitInfo );

            % Create Optimization Problem
            fitProblem = createOptimProblem( 'lsqnonlin', ...
                'objective', f, ...
                'x0', fitVecs.vec, ...
                'ub', fitVecs.ub, ...
                'lb', fitVecs.lb, ...
                'options', opts);

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
                opts = optimoptions( opts, 'OutputFcn', {@plotFit, @fillParams} );
            else
                opts = optimoptions( opts, 'OutputFcn', @fillParams );
            end

            % Obtain the features to fit, their fit vector and labels 
            % Get bounds of parameters to restrict the parameter space
            [ fitVec, fitLabels, ub, lb] = getVec( obj.feature);
            fitVecs = FitEngine.GetVectorBounds(obj, fitVec,fitLabels, ub, lb);

            fitObj = obj.feature; cFeature = 1; 
            fitInfo.Nnew = obj.feature.getSubFeatureNumber();            

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
            fitInfo.featureCurrent.fillParams( size( fitInfo.image));
            fitInfo.timeReversal = obj.parameters.timeReversal;

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
            function stop = fillParams( x, optimV, state)
                FillParams(fitInfo);
                stop = false;
                function FillParams(fitInfo)
                    fitInfo.featureCurrent.fillParams( size( fitInfo.image));
                end
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
                opts = optimoptions( opts, 'OutputFcn', {@plotFit, @fillParams} );
            else
                opts = optimoptions( opts, 'OutputFcn', @fillParams );
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
            function stop = fillParams( x, optimV, state)
                FillParams(fitInfo);
                stop = false;
                function FillParams(fitInfo)
                    fitInfo.featureCurrent.fillParams( size( fitInfo.image));
                end
            end

        end
        % }}}
        
        % }}}
        
        % updateFinalFeature {{{
        function obj = updateFinalFeature( obj, fitInfo)
            
            labels = fitInfo.fitVecs.labels;
            vec = fitInfo.fitResults.vfit;
            if obj.parameters.fitExploreSpeed
                vec = vec .* fitInfo.speedVec;
            end
            fitInfo.featureCurrent.absorbVec(vec, labels);
            
        end
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
                                'MaxIter', 10, ...
                                'TolFun', 1e-9, ...
                                'FiniteDifferenceStepSize', 1e-3, ...
                                'FiniteDifferenceType', 'central', ...
                                'StepTolerance', 1e-5);

            % Set configurable options
            switch config.state
                case 'RELEASE'
                    opts = optimoptions( opts, ...
                                        'display', 'iter',...
                                        'MaxIter', 10);

                case 'DEBUG'
                    opts = optimoptions( opts, ...
                                        'display', 'iter' ,...
                                        'MaxIter', 20);
            end
            
            % Set Parallel Optimization
            if config.useParallel
                opts = optimoptions( opts, 'UseParallel', true);
            end
            
        end
        % }}}

        % SimulateImage{{{
        function [imageOut,err] = SimulateImage( vec, fitInfo)
            % Simulates an image with the feature from fitInfo
            
            % Apply the speedExploration vector to transform to real parameters
            if isfield( fitInfo, 'speedVec')
                try
                vec = vec .* fitInfo.speedVec;
                catch
                    stoph = 1;
                end
            end

            % Get feature to simulate
            feature = fitInfo.featureCurrent;
            featureID = feature.ID;

            % Update feature with vector parameters
            feature.absorbVec( vec, fitInfo.fitVecs.labels );
            
            % Smulate the feature
            % Check if optimization of environment only
            if isfield(fitInfo, 'optimizeEnvOnly') && fitInfo.optimizeEnvOnly == 1
                [imageOut,err] = fitInfo.featureMain.simulateAll( 0*fitInfo.image, 'env');
            else
                [imageOut,err] = fitInfo.featureMain.simulateAll( 0*fitInfo.image, featureID);
            end
            
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
                [imSim, errScale] = FitEngine.SimulateImage( p, fitInfo);
                err = sqrt(errScale) *( imSim - ImageFit);
                err = err( fitInfo.mask(:) );
%                 err = err(:);
            end

        end
        % }}}
        
        % GetVectorBounds {{{
        function fitVecs = GetVectorBounds( obj, vec, vecLabels, varargin)
            % Get upper and lower bounds for the parameters in the fit vector
            
            if length(varargin) == 0
                ub = vec; lb = vec; bs = 1;
            elseif length( varargin) == 2
                ub = varargin{1}; lb = varargin{2}; bs = 0;
            else
                error('unknown number of input arguments')
            end
            
            image = obj.GetImage();
            dim = length( size(image) );

            % Spatial
            minVox = 1;
            maxVox = size( image);
            VoxSize = [maxVox(2), maxVox(1)];
            if dim == 3
                VoxSize = [VoxSize, maxVox(3)];
            end

            % Background
            estBkg = median( image( image > 0) );
            maxBkg = max( image(:) );
            minBkg = min( image(:) );

            % Amplitude
            maxAmp = max( image(:) );
            minAmp = (maxBkg-estBkg )/5; 

            % Sigma
            minSig = [1.2, 1.2, 1.0];
            maxSig = [1.6, 1.6, 1.5];
            if dim == 2
                minSig(end) = []; maxSig(end) = [];
            end

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
                lb( idxAmp) = minAmp; 
                vec( idxAmp( vec(idxAmp) < minAmp) ) = minAmp; end
%             if ~isempty( idxSig), 
%                 nF = length( idxSig)/dim;
%                 ub( idxSig) = repmat( maxSig, 1, nF);
%                 lb( idxSig) = repmat( minSig, 1, nF); end
            if ~isempty( idxP0), 
                nF = length( idxP0)/dim;
                ub( idxP0 ) = repmat( VoxSize, 1, nF);
                lb( idxP0 ) = minVox; end
            if ~isempty( idxP1), 
                nF = length( idxP1)/dim;
                ub( idxP1 ) = repmat( VoxSize, 1, nF);
                lb( idxP1 ) = minVox; end
            if ~isempty( idxP), 
                nF = length( idxP)/dim;
                ub( idxP ) = repmat( VoxSize, 1, nF);
                lb( idxP ) = minVox; end
            if ~isempty( idxBkg), 
                ub( idxBkg ) = maxBkg;
                lb( idxBkg ) = minBkg; end

            % Bounds for Special Objects
            % Interphase curves
            if bs && strcmp(obj.feature.type, 'IMTBank')
                %[ub, lb] = FitEngine.GetVectorBoundsCurves( obj, vec, ub,lb, vecLabels);
                [ub, lb] = FitEngine.GetVectorBoundsBundlesNew( obj, vec, ub,lb, vecLabels);
            end
            if bs && strcmp(obj.feature.type, 'MonopolarAster')
                [ub, lb] = FitEngine.GetVectorBoundsMonopolar( obj, vec, ub,lb, vecLabels);
            end
            
            if any( lb > ub) || any(ub < lb) || any( vec < lb) || any(vec > ub)
                badIdx = unique([ find( lb > ub),  find(ub < lb) , find( vec < lb) , find(vec > ub)]);
                disp('Bad prop labels')
                badProps = vecLabels( badIdx)
                disp(vecLabels( badIdx))
                disp( ub( badIdx) )
                disp( vec( badIdx))
                disp( lb( badIdx) )
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
            % The smaller, the faster
            % Exlporation Speed : unassigned speeds are kept at 1.0
%             speedAmp = 100;
%             speedBkg = 10;
%             speedSigma = 1;
%             speedPos = 10;
            speedAmp = 1;
            speedBkg = 1000;
            speedSigma = 10;
            speedPos = 10;
            speedCXYZ = 10;
            speedT = 10;
            speedL = 100;
            speedNormal = 0.01;
            speedEF = 100;
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
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'cX') ) );
            speedVec( idxAmp) = speedCXYZ;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'cY') ) );
            speedVec( idxAmp) = speedCXYZ;
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'cZ') ) );
            speedVec( idxAmp) = speedCXYZ;
            
            % T
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'origin') ) );speedVec( idxAmp) = speedPos;
            idxT = find( ~cellfun( @isempty, strfind( vecLabels, 'T') ) );speedVec( idxT) = speedT;
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'L') ) );speedVec( idx) = speedT;
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'normalVec') ) ); speedVec( idx) = speedNormal;
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'thetaInit') ) ); speedVec( idx) = speedNormal;
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'ef') ) ); speedVec( idx) = speedEF;
            
            % Monopolar Aster
            % Length
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'length') ) );speedVec( idx) = speedL;
            % Theta
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'theta') ) );speedVec( idx) = 10;
            % Phi
            idx = find( ~cellfun( @isempty, strfind( vecLabels, 'phi') ) );speedVec( idx) = 10;
            
        end
        % }}}
        
    end

end
