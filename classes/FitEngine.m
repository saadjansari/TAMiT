% This class contains the fit machinery required to fit features to an image 
classdef FitEngine

    properties (Access = private)
        image
        image_uint
        feature
        parameters
    end

    methods (Access = public)
        
        % FitEngine {{{
        function obj = FitEngine( image, feature, parameters)
            % Initialize fit Engine with parameters

            obj.image_uint = image;
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
            
            pad_xy = 5;
            pad_z = 1;
            
            % Environment Optimization
            fitInfo = obj.OptimizeEnvironment([]); fitInfo = [];
            
            % Add z-padding
            [obj,~] = obj.zpadding_on( [], pad_xy, pad_z);
            
            fitInfo = [];
            if obj.parameters.fit2DFirst == 1
                [obj,fitInfo] = obj.OptimizeProjection2D();
            end            
            
            % Local Optimization
            fitInfo = obj.OptimizeLocal(fitInfo);
            
            % Global Optimization
            fitInfo = obj.OptimizeGlobal(fitInfo);
            
            % Hyper Parameters Optimization
            fitInfo = obj.OptimizeHyperParameters( fitInfo);
            
            % Save final optimized data
            fitInfo.fitScope = 'final';
            obj.feature = fitInfo.featureMain;
            
            % initial/final images
            imageSimI = FitEngine.SimulateImage( fitInfo.fitVecs.vec, fitInfo);
            imageSimF = FitEngine.SimulateImage( fitInfo.fitResults.vfit, fitInfo);
            fitInfo.imageSimI = imageSimI;
            fitInfo.imageSimF = imageSimF;
            
            % Remove z-padding
            [obj,fitInfo] = obj.zpadding_off( fitInfo, pad_xy, pad_z);
            
            % Scale amplitudes
%             [obj, fitInfo] = obj.scale_amplitude( fitInfo);
%             fitInfo.featureMain.image = obj.image_uint;
            
            obj.feature.finalizeAddedFeatures();
            if obj.parameters.display
                Cell.displayFinalFit( obj.image, obj.feature, fitInfo);
            end
            Cell.saveFinalFit( obj.image_uint, obj.feature, fitInfo);


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
                
                try
                    % Prepare for optimization
                    obj.feature.forceInsideMask();
                    [ fitProblem{jFeature}, fitInfo{jFeature}] = obj.PrepareOptimizeLocal( jFeature);

                    % Solve Optimization Problem
                    fitInfo{ jFeature}.fitResults = obj.SolveOptimizationProblem( fitProblem{ jFeature} );
                    obj.updateFinalFeature( fitInfo{jFeature} );

                    % Update feature heirarchy to ensure all dependencies are acknowledged
                    obj.feature.updateSubFeatures();
                catch
                    disp('Could not locally fit feature')
                end

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
            binit = obj.feature.background;
            if ~isempty(obj.feature.backgroundNuclear)
                bninit = obj.feature.backgroundNuclear;
                fprintf('\t Initial = %.4f\n\t Initial Nuclear = %.4f\n',binit, bninit)
                
                blist = linspace(0.7*binit, 1.5*binit, 20);
                bnlist = linspace(0.7*bninit, 1.5*bninit, 20);
                res = zeros( [length(blist), length(bnlist)]);
                for ib = 1 : length(blist)
                    for ibn = 1 : length(bnlist)
                        obj.feature.background = blist(ib);
                        obj.feature.backgroundNuclear = bnlist(ibn);
                        iSim = obj.feature.simulateAll( 0*obj.feature.image, obj.feature.ID);
                        res(ib,ibn) = sum( (obj.feature.image(:) - iSim(:)).^2 );
                    end
                end
                [~, idx0] = min( min( res,[],2) );
                [~, idx1] = min( min( res,[],1) );
                obj.feature.background = blist(idx0);
                obj.feature.backgroundNuclear = bnlist(idx1);
                fprintf('\t Final = %.4f\n\t Final Nuclear = %.4f\n',obj.feature.background, obj.feature.backgroundNuclear)
            else
                fprintf('  Initial = %.4f\n',binit)
                blist = linspace(0.5*binit, 2*binit, 100);
                res = zeros( size(blist));
                for ib = 1 : length(blist)
                    obj.feature.background = blist(ib);
                    iSim = obj.feature.simulateAll( 0*obj.feature.image, obj.feature.ID);
                    res(ib) = sum( (obj.feature.image(:) - iSim(:)).^2 );
                end
                [~, idx] = min( res);
                obj.feature.background = blist(idx);
                fprintf('  Final = %.4f\n', obj.feature.background)
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
                refImage = obj.image - FitEngine.SimulateImage( fitInfoOld.fitResults.vfit, fitInfoOld);
               
                % Create a deep copy main feature
                featureNew = feature.copyDeep();

                % Add Features
                switch featureType
                    case 'Basic'
                        [~, successAdd ] = featureNew.addSubFeatures( refImage);
                    case 'Organizer'
                        [~, successAdd ] = featureNew.addOrganizers( refImage);
                end
                
                featureNew.syncFeaturesWithMap();
                
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
                    
                    % Update
                    labels = fitInfo.fitVecs.labels;
                    vec = fitInfo.fitResults.vfit;
                    vecErr = fitInfo.fitResults.vfitError;
                    if obj.parameters.fitExploreSpeed
                        vec = vec .* fitInfo.speedVec;
                        vecErr = vecErr .* fitInfo.speedVec;
                    end
                    featureNew.absorbVec(vec, labels);
                    featureNew.absorbVec(vecErr, labels, 1);

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
                featureNew.syncFeaturesWithMap();
                
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
                    
                    % Update
                    labels = fitInfo.fitVecs.labels;
                    vec = fitInfo.fitResults.vfit;
                    vecErr = fitInfo.fitResults.vfitError;
                    if obj.parameters.fitExploreSpeed
                        vec = vec .* fitInfo.speedVec;
                        vecErr = vecErr .* fitInfo.speedVec;
                    end
                    featureNew.absorbVec(vec, labels);
                    featureNew.absorbVec(vecErr, labels, 1);
                    
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
            % Use to calculate confidence intervals
            % This replaces nlparci 
            estParams = fitResults.vfit;
            residual = fitResults.residual;
            jacobian = fitResults.jacobian;
            A = full(jacobian);
            B = nearestSPD(A'*A);
            % fprintf('Determinant = %.3e\n',det(B))

            lastwarn(''); %Clear warning memory
            alpha = 0.05; %significance level
            df = length(residual) - numel(estParams); %degrees of freedom
            crit = tinv(1-alpha/2,df);       %critical value
            % covm = inv(jacobian'*jacobian) * var(residual); %covariance matrix
            covm = B; %covariance matrix
            [~, warnId] = lastwarn; %detect if inv threw 'nearly singular' warning.
            covmIdx = sub2ind(size(covm),1:size(covm,1),1:size(covm,2));  %indices of the diag of covm
            CI = nan(numel(estParams),2);
            if ~strcmp(warnId, 'MATLAB:nearlySingularMatrix')
                CI(:,1) = estParams - crit * sqrt(covm(covmIdx));
                CI(:,2) = estParams + crit * sqrt(covm(covmIdx));
            end
            fitResults.vfitError = abs( (CI(:,2)-CI(:,1) ) / 2 )';

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
            
            % Update 3D features using fitted 2D information
            obj.feature.Update3DFrom2D(fitInfo.featureCurrent);
            obj.feature.syncFeatures();
            
            obj.feature.fit = 'zonly';
            for j1 = 1 : obj.feature.numFeatures
                obj.feature.featureList{j1}.fit = 'zonly';
                if isprop( obj.feature.featureList{j1}, 'featureList')
                    for j2 = 1 : obj.feature.featureList{j1}.numFeatures
                        obj.feature.featureList{j1}.featureList{j2}.fit = 'zonly';
                        if isprop( obj.feature.featureList{j1}.featureList{j2}, 'featureList')
                            for j3 = 1 : obj.feature.featureList{j1}.featureList{j2}.numFeatures
                                obj.feature.featureList{j1}.featureList{j2}.featureList{j3}.fit = 'zonly';
                            end
                        end
                    end
                end
            end
            
            % Global Z-Optimization
            fitInfo = obj.OptimizeGlobal(fitInfo);
            
            if isfield( fitInfo, 'speedVec')
                vec = vec .* fitInfo.speedVec;
                fitInfo.featureCurrent.absorbVec( fitInfo.fitResults.vfit.*fitInfo.speedVec, fitInfo.fitVecs.labels);
            end
            
            % Unscale parameters
            if isfield( fitInfo.fitVecs, 'scaleParameters')
                if fitInfo.fitVecs.scaleParameters == 1
                    vec_unscaled = FitEngine.unscale_parameters( fitInfo.fitResults.vfit, fitInfo.fitVecs.ub, fitInfo.fitVecs.lb, fitInfo.fitVecs.ub_unscaled, fitInfo.fitVecs.lb_unscaled);
                    fitInfo.featureCurrent.absorbVec( vec_unscaled, fitInfo.fitVecs.labels);
                end
            end
            
            obj.feature.fit = 'all';
            for j1 = 1 : obj.feature.numFeatures
                obj.feature.featureList{j1}.fit = 'all';
                if isprop( obj.feature.featureList{j1}, 'featureList')
                    for j2 = 1 : obj.feature.featureList{j1}.numFeatures
                        obj.feature.featureList{j1}.featureList{j2}.fit = 'all';
                        if isprop( obj.feature.featureList{j1}.featureList{j2}, 'featureList')
                            for j3 = 1 : obj.feature.featureList{j1}.featureList{j2}.numFeatures
                                obj.feature.featureList{j1}.featureList{j2}.featureList{j3}.fit = 'all';
                            end
                        end
                    end
                end
            end
            
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
                speedVec = obj.getExplorationSpeedVector( fitLabels{ jFeature}, obj.feature.type);
                fitVecs.vec = fitVecs.vec ./ speedVec;
                fitVecs.ub = fitVecs.ub ./ speedVec;
                fitVecs.lb = fitVecs.lb ./ speedVec;
                fitInfo.speedVec = speedVec;
            end
            if obj.parameters.scaleParameters
                fitVecs.scaleParameters = 1;
                [vec_scaled, ub_scaled, lb_scaled] = FitEngine.scale_parameters( fitVecs.vec, fitVecs.ub, fitVecs.lb);
                fitVecs.vec_unscaled = fitVecs.vec;
                fitVecs.ub_unscaled = fitVecs.ub;
                fitVecs.lb_unscaled = fitVecs.lb;
                fitVecs.vec = vec_scaled;
                fitVecs.ub = ub_scaled;
                fitVecs.lb = lb_scaled;
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
                speedVec = obj.getExplorationSpeedVector( fitLabels, obj.feature.type);
                fitVecs.vec = fitVecs.vec ./ speedVec;
                fitVecs.ub = fitVecs.ub ./ speedVec;
                fitVecs.lb = fitVecs.lb ./ speedVec;
                fitInfo.speedVec = speedVec;
            end
            if obj.parameters.scaleParameters
                fitVecs.scaleParameters = 1;
                [vec_scaled, ub_scaled, lb_scaled] = FitEngine.scale_parameters( fitVecs.vec, fitVecs.ub, fitVecs.lb);
                fitVecs.vec_unscaled = fitVecs.vec;
                fitVecs.ub_unscaled = fitVecs.ub;
                fitVecs.lb_unscaled = fitVecs.lb;
                fitVecs.vec = vec_scaled;
                fitVecs.ub = ub_scaled;
                fitVecs.lb = lb_scaled;
            end
            
            % Create a parameter tracker
            p_tracker = fitVecs.vec;
            
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
            fitInfo.param_tracker = p_tracker;

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
                % track parameters
                fitInfo.param_tracker = [fitInfo.param_tracker ; x];
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
            vecErr = fitInfo.fitResults.vfitError;
            if obj.parameters.fitExploreSpeed
                vec = vec .* fitInfo.speedVec;
                vecErr = vecErr .* fitInfo.speedVec;
            end
            fitInfo.featureCurrent.absorbVec(vec, labels);
            fitInfo.featureCurrent.absorbVec(vecErr, labels, 1);
        end
        % }}}
        
        function [obj,fitInfo] = zpadding_on( obj, fitInfo, pad_xy, pad_z)
           
            % Pad image
            img0 = obj.image;
            img1 = zeros( 2*pad_xy+size(img0, 1), 2*pad_xy+size(img0,2), 2*pad_z+size(img0,3), class(img0) );
            img1( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z) = img0;
            obj.image = img1;
            
            % Pad features
            % Main organizer
            obj.feature.image = img1;
            mask0 = obj.feature.mask;
            mask1 = zeros( size(img1,1), size(img1,2), size(img1,3), class(mask0) );
            mask1( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z) = mask0;
            mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1) = mask0(:,:,1);
            mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, end) = mask0(:,:,end);
            obj.feature.mask = mask1;
            
            if ~isempty(obj.feature.maskNuclear)
                mask0 = obj.feature.maskNuclear;
                mask1 = zeros( size(img1,1), size(img1,2), size(img1,3), class(mask0) );
                mask1( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z) = mask0;
                mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1) = mask0(:,:,1);
                mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, end) = mask0(:,:,end);
                obj.feature.maskNuclear = mask1;
            end
            
            switch obj.feature.type
                case 'MonopolarAster'
                    faster = obj.feature.featureList{1};
                    for jf = 1 : faster.numFeatures
                        
                        if strcmp( faster.featureList{jf}.type, 'Spot')
%                             faster.featureList{jf}.position(3) = faster.featureList{jf}.position(3)+1;
%                             faster.featureList{jf}.bounds.ub.position(3) = faster.featureList{jf}.bounds.ub.position(3) + 2;
                            faster.featureList{jf}.position = faster.featureList{jf}.position + [pad_xy,pad_xy,pad_z];
                            faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position + 2*[pad_xy,pad_xy,pad_z];
                        elseif strcmp( faster.featureList{jf}.type, 'Line')
%                             faster.featureList{jf}.startPosition(3) = faster.featureList{jf}.startPosition(3)+1;
%                             faster.featureList{jf}.endPosition(3) = faster.featureList{jf}.endPosition(3)+1;
%                             faster.featureList{jf}.bounds.ub.startPosition(3) = faster.featureList{jf}.bounds.ub.startPosition(3) + 2;
%                             faster.featureList{jf}.bounds.ub.endPosition(3) = faster.featureList{jf}.bounds.ub.endPosition(3) + 2;
                            faster.featureList{jf}.startPosition = faster.featureList{jf}.startPosition + [pad_xy,pad_xy,pad_z];
                            faster.featureList{jf}.endPosition = faster.featureList{jf}.endPosition + [pad_xy,pad_xy,pad_z];
                            faster.featureList{jf}.bounds.ub.startPosition = faster.featureList{jf}.bounds.ub.startPosition + 2*[pad_xy,pad_xy,pad_z];
                            faster.featureList{jf}.bounds.ub.endPosition = faster.featureList{jf}.bounds.ub.endPosition + 2*[pad_xy,pad_xy,pad_z];
                        else
                            error('unknown feature type')
                        end
                        faster.featureList{jf}.fillParams( size(img1) );
                    end
                    
                case 'Spindle'
                    
                    fspindle = obj.feature.featureList{1};
                    fspindle.startPosition = fspindle.startPosition + [pad_xy,pad_xy,pad_z];
                    fspindle.endPosition = fspindle.endPosition + [pad_xy,pad_xy,pad_z];
                    fspindle.bounds.ub.startPosition = fspindle.bounds.ub.startPosition + 2*[pad_xy,pad_xy,pad_z];
                    fspindle.bounds.ub.endPosition = fspindle.bounds.ub.endPosition + 2*[pad_xy,pad_xy,pad_z];
                    
                    for ja = 2 : length( obj.feature.featureList)
                        faster =  obj.feature.featureList{ja};
                        for jf = 1 : faster.numFeatures

                            if strcmp( faster.featureList{jf}.type, 'Spot')
                                faster.featureList{jf}.position = faster.featureList{jf}.position + [pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position + 2*[pad_xy,pad_xy,pad_z];
                            elseif strcmp( faster.featureList{jf}.type, 'Line')
    %                             faster.featureList{jf}.startPosition(3) = faster.featureList{jf}.startPosition(3)+1;
    %                             faster.featureList{jf}.endPosition(3) = faster.featureList{jf}.endPosition(3)+1;
    %                             faster.featureList{jf}.bounds.ub.startPosition(3) = faster.featureList{jf}.bounds.ub.startPosition(3) + 2;
    %                             faster.featureList{jf}.bounds.ub.endPosition(3) = faster.featureList{jf}.bounds.ub.endPosition(3) + 2;
                                faster.featureList{jf}.startPosition = faster.featureList{jf}.startPosition + [pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.endPosition = faster.featureList{jf}.endPosition + [pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.bounds.ub.startPosition = faster.featureList{jf}.bounds.ub.startPosition + 2*[pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.bounds.ub.endPosition = faster.featureList{jf}.bounds.ub.endPosition + 2*[pad_xy,pad_xy,pad_z];
                            else
                                error('unknown feature type')
                            end
                            faster.featureList{jf}.fillParams( size(img1) );
                        end
                    end
                    
                case 'SpindleNew'
                    
                    fspindle = obj.feature.featureList{1};
                    fspindle.startPosition = fspindle.startPosition + [pad_xy,pad_xy,pad_z];
                    fspindle.endPosition = fspindle.endPosition + [pad_xy,pad_xy,pad_z];
                    fspindle.bounds.ub.startPosition = fspindle.bounds.ub.startPosition + 2*[pad_xy,pad_xy,pad_z];
                    fspindle.bounds.ub.endPosition = fspindle.bounds.ub.endPosition + 2*[pad_xy,pad_xy,pad_z];
                    
                    for ja = 2 : length( obj.feature.featureList)
                        faster =  obj.feature.featureList{ja};
                        for jf = 1 : faster.numFeatures

                            if strcmp( faster.featureList{jf}.type, 'Spot')
                                faster.featureList{jf}.position = faster.featureList{jf}.position + [pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position + 2*[pad_xy,pad_xy,pad_z];
                            elseif strcmp( faster.featureList{jf}.type, 'CurvedMT')
                                faster.featureList{jf}.origin = faster.featureList{jf}.origin + [pad_xy,pad_xy,pad_z];
                                faster.featureList{jf}.bounds.ub.origin = faster.featureList{jf}.bounds.ub.origin + 2*[pad_xy,pad_xy,pad_z];
                            else
                                error('unknown feature type')
                            end
                            faster.featureList{jf}.fillParams( size(img1) );
                        end
                    end
                    
            end
        end
        
        function [obj,fitInfo] = zpadding_off( obj, fitInfo, pad_xy, pad_z)
           
            % unPad image
            obj.image = obj.image( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            
            % unPad features
            obj.feature = fitInfo.featureMain;
            % Main organizer
            obj.feature.image = obj.feature.image( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            obj.feature.mask = obj.feature.mask( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            if ~isempty(obj.feature.maskNuclear)
                obj.feature.maskNuclear = obj.feature.maskNuclear( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            end
            obj.feature.imageSim = obj.feature.imageSim( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            
            %fitInfo
            fitInfo.imageSimI = fitInfo.imageSimI( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            fitInfo.imageSimF = fitInfo.imageSimF( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
            
            pad_arr = [pad_xy, pad_xy, pad_z];
            switch obj.feature.type
                case 'MonopolarAster'
                    faster = obj.feature.featureList{1};
                    %faster.imageSim = faster.imageSim( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z);
                    
                    for jf = 1 : faster.numFeatures
                        if strcmp( faster.featureList{jf}.type, 'Spot')
                            faster.featureList{jf}.position = faster.featureList{jf}.position - pad_arr;
                            faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position - 2*pad_arr;
                            
                        elseif strcmp( faster.featureList{jf}.type, 'Line')
                            faster.featureList{jf}.startPosition = faster.featureList{jf}.startPosition - pad_arr;
                            faster.featureList{jf}.endPosition = faster.featureList{jf}.endPosition - pad_arr;
                            faster.featureList{jf}.bounds.ub.startPosition = faster.featureList{jf}.bounds.ub.startPosition - 2*pad_arr;
                            faster.featureList{jf}.bounds.ub.endPosition = faster.featureList{jf}.bounds.ub.endPosition - 2*pad_arr;
                        else
                            error('unknown feature type')
                        end
                        faster.featureList{jf}.fillParams( size(obj.image) );
                    end
                    
                case 'Spindle'
                    
                    fspindle = obj.feature.featureList{1};
                    fspindle.startPosition = fspindle.startPosition - pad_arr;
                    fspindle.endPosition = fspindle.endPosition - pad_arr;
                    fspindle.bounds.ub.startPosition = fspindle.bounds.ub.startPosition - 2*pad_arr;
                    fspindle.bounds.ub.endPosition = fspindle.bounds.ub.endPosition - 2*pad_arr;
                    
                    for ja = 2 : length( obj.feature.featureList)
                        faster =  obj.feature.featureList{ja};
                        for jf = 1 : faster.numFeatures

                            if strcmp( faster.featureList{jf}.type, 'Spot')
                                faster.featureList{jf}.position = faster.featureList{jf}.position - pad_arr;
                                    faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position  - 2*pad_arr;
                            elseif strcmp( faster.featureList{jf}.type, 'Line')
                                faster.featureList{jf}.startPosition = faster.featureList{jf}.startPosition - pad_arr;
                                faster.featureList{jf}.endPosition = faster.featureList{jf}.endPosition - pad_arr;
                                faster.featureList{jf}.bounds.ub.startPosition = faster.featureList{jf}.bounds.ub.startPosition - 2*pad_arr;
                                faster.featureList{jf}.bounds.ub.endPosition = faster.featureList{jf}.bounds.ub.endPosition - 2*pad_arr;
                            else
                                error('unknown feature type')
                            end
                            faster.featureList{jf}.fillParams( size(obj.image) );
                        end
                    end
                    
                    
                case 'SpindleNew'
                    
                    fspindle = obj.feature.featureList{1};
                    fspindle.startPosition = fspindle.startPosition - pad_arr;
                    fspindle.endPosition = fspindle.endPosition - pad_arr;
                    fspindle.bounds.ub.startPosition = fspindle.bounds.ub.startPosition - 2*pad_arr;
                    fspindle.bounds.ub.endPosition = fspindle.bounds.ub.endPosition - 2*pad_arr;
                    
                    for ja = 2 : length( obj.feature.featureList)
                        faster =  obj.feature.featureList{ja};
                        for jf = 1 : faster.numFeatures

                            if strcmp( faster.featureList{jf}.type, 'Spot')
                                faster.featureList{jf}.position = faster.featureList{jf}.position - pad_arr;
                                faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position  - 2*pad_arr;
                            elseif strcmp( faster.featureList{jf}.type, 'CurvedMT')
                                faster.featureList{jf}.origin = faster.featureList{jf}.origin  - pad_arr;
                                faster.featureList{jf}.bounds.ub.origin = faster.featureList{jf}.bounds.ub.origin  - 2*pad_arr;
                            else
                                error('unknown feature type')
                            end
                            faster.featureList{jf}.fillParams( size(obj.image) );
                        end
                    end
            end
            
            % Alter fitInfo
            fitInfo.featureMain = obj.feature;
            fitInfo.mask = obj.feature.mask;
            fitInfo.image = obj.feature.image;
            fitInfo.numVoxels = size(fitInfo.image);

        end
        
        function [obj, fitInfo] = scale_amplitude( obj, fitInfo)
            
            imraw = obj.image_uint;
            img = fitInfo.featureMain.image;
            
            scale_factor = max(imraw(:)) / max(img(:));
            
            switch fitInfo.featureMain.type
                case 'MonopolarAster'
                    f1 = fitInfo.featureMain.featureList{1};
                    for j1 = 1 : length( f1.featureList)
                        f1.featureList{j1}.amplitude = scale_factor * f1.featureList{j1}.amplitude;
                        f1.featureList{j1}.err_amplitude = scale_factor * f1.featureList{j1}.err_amplitude;
                    end                    
                case 'Mitosis'
                    error('not set up')
                    f1 = fitInfo.featureMain.featureList{1};
                    f2 = fitInfo.featureMain.featureList{2};
                    f3 = fitInfo.featureMain.featureList{3};
                    for j1 = 1 : length( f1.featureList)
                        f1.featureList{j1}.amplitude = scale_factor * f1.featureList{j1}.amplitude;
                        f1.featureList{j1}.err_amplitude = scale_factor * f1.featureList{j1}.err_amplitude;
                    end  
                    
                case 'MitosisBud'
                    
            end
            
        end
        
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
                                'FiniteDifferenceStepSize', 1e-5, ...
                                'FiniteDifferenceType', 'central', ...
                                'StepTolerance', 1e-5);

            % Set configurable options
            switch config.state
                case 'RELEASE'
                    opts = optimoptions( opts, ...
                                        'display', 'iter',...
                                        'MaxIter', 30);

                case 'DEBUG'
                    opts = optimoptions( opts, ...
                                        'display', 'iter' ,...
                                        'MaxIter', 30);
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
            
            % Unscale parameters
            if isfield( fitInfo.fitVecs, 'scaleParameters')
                if fitInfo.fitVecs.scaleParameters == 1
                    vec_unscaled = FitEngine.unscale_parameters( vec, fitInfo.fitVecs.ub, fitInfo.fitVecs.lb, fitInfo.fitVecs.ub_unscaled, fitInfo.fitVecs.lb_unscaled);
                    vec = vec_unscaled;
                end
            end

            % Get feature to simulate
            feature = fitInfo.featureCurrent;
            featureID = feature.ID;

            % Update feature with vector parameters
            feature.absorbVec( vec, fitInfo.fitVecs.labels );
            
            % Smulate the feature
            % Check if optimization of environment only
            [imageOut,err] = fitInfo.featureMain.simulateAll( 0*fitInfo.image, featureID);
            
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
                err = 1 *( imSim - ImageFit);
%                 err = err( fitInfo.mask(:) );
                err = err(:);
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
            minBkg = min( [min( image(:) ), 0]);

            % Find the correct label in vecLabels, and place the correct bounds in the correct places

            % Find index of parameters 
            idxP0 = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            idxP1 = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            idxP = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            idxOrigin = find( ~cellfun( @isempty, strfind( vecLabels, 'origin') ) );
            idxBkg = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );
            
            % Store upper and lower bounds correctly
            if strcmp( obj.feature.fit, 'zonly')
                VoxSize = VoxSize(end); dim=1;
            end
            
            %positions
            if ~isempty(idxP0)
                nF = length(idxP0)/dim;
                ub(idxP0) = repmat( VoxSize, 1, nF);
                lb(idxP0) = minVox;
            end
            if ~isempty(idxP1)
                nF = length(idxP1)/dim;
                ub(idxP1) = repmat( VoxSize, 1, nF);
                lb(idxP1) = minVox;
            end
            if ~isempty(idxP)
                nF = length(idxP)/dim;
                ub(idxP) = repmat( VoxSize, 1, nF);
                lb(idxP) = minVox;
            end
            if ~isempty(idxOrigin)
                nF = length(idxOrigin)/dim;
                ub(idxOrigin) = repmat( VoxSize, 1, nF);
                lb(idxOrigin) = minVox;
            end
            if ~isempty(idxBkg)
                ub(idxBkg) = maxBkg;
                lb(idxBkg) = minBkg;
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
        function scaleVec = getExplorationSpeedVector( vecLabels, featType)
            % Creates a weighing vector to allow a user to assign different importance to different kinds of parameters. This will infact allow the fitting optimization engine to explore the parameters space at different speeds
            % The parameters are divided by this vector
            % Exlporation Speed : unassigned speeds are kept at 1.0
            
            % NOTE NOTE NOTE!!!
            % This should probably be output directly from the main
            % feature.
            
            scaleVec = ones( size(vecLabels) );

            % function to find indices of matching substrings
            find_str = @(x) find( ~cellfun( @isempty, strfind( vecLabels, x) ) );
            
            % Common Scalings
            scale.amplitude = 0.01;
            scale.background = 0.1;
            scale.sigma = 0.01;
            scale.position = 0.1;
            
            scaleVec( find_str('amplitude') ) = scale.amplitude;
            scaleVec( find_str('background') ) = scale.background;
            scaleVec( find_str('sigma') ) = scale.sigma;
            scaleVec( find_str('position') ) = scale.position;
            scaleVec( find_str('startPosition') ) = scale.position;
            scaleVec( find_str('endPosition') ) = scale.position;
            
            
            
            switch featType
                case 'MonopolarAster'
                    scale.theta = 0.1;
                    scale.length = 1;
                    scaleVec( find_str('theta') ) = scale.theta;
                    scaleVec( find_str('length') ) = scale.length;
                    
                case 'Spindle'
                    scale.theta = 0.1;
                    scale.length = 1;
                    scaleVec( find_str('theta') ) = scale.theta;
                    scaleVec( find_str('length') ) = scale.length;
                    
                case 'SpindleNew'
                    scale.normal_vec = 0.001;
                    scale.theta = 0.01;
                    scale.length = 10;
                    scaleVec( find_str('thetaInit') ) = scale.theta;
                    scaleVec( find_str('L') ) = scale.length;
                    scaleVec( find_str('normalVec') ) = scale.normal_vec;
                    scaleVec( find_str('origin') ) = scale.position;
                otherwise
                    error('unknown feature type')
                    
            end

            
            % T
%             speedCXYZ = 10;
%             speedT = 10;
%             speedL = 1;
%             speedNormal = 0.01;
%             speedTheta = 0.01;
%             speedEF = 100;
%             idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, '_origin') ) );scaleVec( idxAmp) = speedPos;
%             idxT = find( ~cellfun( @isempty, strfind( vecLabels, '_T') ) );scaleVec( idxT) = speedT;
%             idxL = find( ~cellfun( @isempty, strfind( vecLabels, '_L') ) );scaleVec( idxL) = speedL;
%             idxN = find( ~cellfun( @isempty, strfind( vecLabels, 'normalVec') ) ); scaleVec( idxN) = speedNormal;
%             idxTheta = find( ~cellfun( @isempty, strfind( vecLabels, 'thetaInit') ) ); scaleVec( idxTheta) = speedTheta;
%             idxEF = find( ~cellfun( @isempty, strfind( vecLabels, '_ef') ) ); scaleVec( idxEF) = speedEF;
            
%             % Monopolar Aster
%             % Length
%             idx = find( ~cellfun( @isempty, strfind( vecLabels, '_length') ) );scaleVec( idx) = speedL;
%             % Theta
%             idx = find( ~cellfun( @isempty, strfind( vecLabels, '_theta') ) );scaleVec( idx) = 10;
%             % Phi
%             idx = find( ~cellfun( @isempty, strfind( vecLabels, '_phi') ) );scaleVec( idx) = 10;
            %speedVec = ones( size(speedVec));
        end
        % }}}
        
        % Scale parameters for fitting {{{
        function [vec_scaled, ub_scaled, lb_scaled] = scale_parameters( vec, ub, lb)
            
            % Ensure vec is between lb and ub
            assert( all(vec >= lb), 'Some initial parameters are below the lower bound.')
            assert( all(vec <= ub), 'Some initial parameters are above the upper bound.')
            assert( all(ub >= lb), 'Some upper bound values are below the lower bound.')
            
            % Lower bound will be 0, Upper bound will be 1.
            lb_scaled = zeros( size(vec) );
            ub_scaled = ones( size(vec) );
            
            % Scale vector
            vec_scaled = (vec-lb)./(ub-lb);
        end
        
        % Unscale parameters for fitting
        function vec_unscaled = unscale_parameters( vec_scaled, ub_scaled, lb_scaled, ub_unscaled, lb_unscaled)
                                    
            % Scale vector
            vec_unscaled = lb_unscaled + (vec_scaled-lb_scaled)./(ub_scaled-lb_scaled) .* (ub_unscaled-lb_unscaled);
            
        end
        % }}}
    end

end
