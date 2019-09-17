classdef Cell < handle & matlab.mixin.Copyable
% This is a generic cell which can be specialized
    properties
        type % type of cell (e.g. interphase or mitotis). THis allows for multiple implementations depending on cell phase
        species
        strain % any specific strain information
        numberOfChannels
        featuresInChannels
        channelsToFit
        lifetime % 2-element vector containing [start_time, end_time]. Should be used when accessing image
        image % (X,Y,Z,T,C)
        featureList % cell array containing features. this is of size numberOfChannels x Time.
        settings
        featureMap
    end

    methods ( Access = public )
        
        % Cell {{{
        function obj = Cell( image, lifetime, species, featuresInChannels, channelsToFit, type, settings)
        % Cell : constructor function for Cell superclass
       
            % Store essentials
            obj.image = Cell.img2double( image);
            obj.lifetime = lifetime;
            obj.species = species;
            obj.featuresInChannels = featuresInChannels;
            obj.channelsToFit = channelsToFit;
            obj.numberOfChannels = size( image, 5);
            obj.settings = Cell.parseSettings( settings);

            % Initialize the feature list
            obj.featureList = cell( length( obj.channelsToFit), size(image, 4) );
            
            % Initialize featureIdxMap
            obj.featureMap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');

            % Cell Type
            if nargin < 6, obj.type = 'generic';
            else, obj.type = type; end

            % Save Directory default behaviour
            if nargin < 7, 
                obj.settings.saveDirectory = [pwd, filesep, 'ResultsDump'];
                warning('off', 'MATLAB:RMDIR:RemovedFromPath')
                rmdir( obj.settings.saveDirectory, 's');
                warning('on', 'MATLAB:RMDIR:RemovedFromPath')
            end
            if exist( obj.settings.saveDirectory) ~= 7
                mkdir( obj.settings.saveDirectory)
            end

            % Displays
            disp('    -----------------------   C E L L    I N F O   ------------------------');
            disp('    -----------------------------------------------------------------------');disp(' ')
            disp( ['        Type : ' obj.type ] )
            disp( ['        Species : ' , obj.species])
            disp( sprintf( '        Lifetime : %d - %d', lifetime(1), lifetime(2) ) );
            disp( sprintf( '        Movie size XYZTC = %d x %d x %d x %d x %d', size(image, 1), size(image, 2), size(image, 3), size(image, 4), size( image, 5) ) )
            disp( ['        Features = ', strjoin( featuresInChannels, ' - ') ] ) 
            disp('     ----------------------------------------------------------------------');disp(' ')

        end
        % }}}

        % fitFeatures {{{
        function obj = fitFeatures( obj, parameters)

            for jChannel = 1 : length(obj.channelsToFit)

                % get actual channel number based on channel index
                cChannel = obj.channelsToFit( jChannel);
                disp( upper( sprintf( ['                               Channel %d = ' obj.featuresInChannels{jChannel}], cChannel ) ) )
                disp('    -----------------------------------------------------------------------')

                for jTime = obj.lifetime(1) : obj.lifetime(2)

                    disp( sprintf( '    C-%d Time = %d', jChannel, jTime) )

                    % Initialize for this CT step
                    parameters.channelTrue = cChannel;
                    parameters.channelIdx = jChannel;
                    parameters.time = jTime;
                    parameters.saveDirectory = [ obj.settings.saveDirectory , filesep, sprintf( 'C%d_T%d', cChannel, jTime) ];
                    mkdir( parameters.saveDirectory);

                    % Fit features for this CT frame
                    obj.fitFeatures_Frame( parameters) 

                end
                disp('    -----------------------------------------------------------------------')
            end
            
        end
        
        % fitFeatures_Frame {{{
        function fitFeatures_Frame( obj, parameters)
            % Find and fit features for a single CT frame
            
            Image2Fit = obj.image(:,:,:, parameters.time, parameters.channelTrue);

            % check if frame should be skipped
            status = Cell.checkFrameViability( Image2Fit, parameters.time);
            if status ~= 0
                obj.featureList{ parameters.channelIdx, parameters.time} = NaN;
                return
            end

            % Estimate the features based on an estimation routine (defined in specialized sub-class )
            % Good estimation is key to good optimization in low SnR images
            disp('        Estimation')
            obj = findFeatures( obj, parameters.time, parameters.channelTrue, parameters.channelIdx); 
            mainFeature = obj.featureList{ parameters.channelIdx , parameters.time};


            % Optimize the features via local and global fitting, followed by optimization of feature number
            disp('        Optimization')
            disp('            - Regular')

            % Local Fitting
            obj.settings.flags.doFitLocal = 0;
            if obj.settings.flags.doFitLocal
                fitInfo.Local = fitLocal(  obj, Image2Fit, parameters);
            end

            % Global Fitting
            fitInfo.Global = fitGlobal( obj, Image2Fit, parameters);

            % Display and Save
            if obj.settings.flags.graphRealTime
                Cell.displayFinalFit( Image2Fit, mainFeature, fitInfo.Global);
            end
            Cell.saveFinalFit( Image2Fit, mainFeature, fitInfo.Global);

            % HyperParameter Fitting
            if obj.settings.flags.doFitFeatureNumber && ~obj.settings.flags.skipOptimizeNumber
                disp('            - Feature Number')
                fitInfo.GlobalNumber = fitFeatureNumber( obj, Image2Fit, parameters, fitInfo.Global);

                % Display and Save
                if obj.settings.flags.graphRealTime
                    Cell.displayFinalFit( Image2Fit, obj.featureList{parameters.channelIdx, parameters.time}, fitInfo.GlobalNumber);
                end
                Cell.saveFinalFit( Image2Fit, mainFeature, fitInfo.GlobalNumber);
            end

            obj.updateFeatureMap( parameters.channelIdx, parameters.time);
            close all
            
            %fitInfo.Global.fitVecs.labels
            %fitInfo.Global.fitVecs.vec
            %fitInfo.Global.fitResults.vfit

        end
        % }}}
        % fitLocal {{{
        function fitInfo = fitLocal( obj, Image2Fit, parameters)

            disp('                    Local Fitting : '); 

            % Check if mainFeature should be fit. If not, then return out of method
            if isempty( mainFeature.featureList)
                fprintf('                       Skipping fitting for %s\n', mainFeature.type);  
                [ ~, fitInfo] = prepareFit( obj, mainFeature, Image2Fit, parameters, 'local', 1);
                fitInfo.featureMain = mainFeature;
                return
            end
            
            mainFeature = obj.featureList{ parameters.channelIdx, parameters.time};
            nFeatures = mainFeature.getSubFeatureNumber();

            for jFeature = 1 : nFeatures
                disp( sprintf('                        Current Feature = %d / %d', jFeature, nFeatures) ); 

                % Prepare for the Fit (specialized function)
                [ fitProblem{jFeature}, fitInfo{jFeature}] = prepareFit( obj, mainFeature, Image2Fit, parameters, 'local', jFeature);

                % run the lsqnonlin local fitting
                [fR.vfit,fR.resnorm,fR.residual,fR.exitflag,~,~,fR.jacobian] = lsqnonlin( fitProblem{ jFeature} ); 
                
                % Get the errors from the jacobian
                ci = nlparci( fR.vfit, fR.residual, 'jacobian', fR.jacobian, 'alpha', 0.07); 
                fR.vfitError = abs( (ci(:,2)-ci(:,1) ) / 2 )';
                fitInfo{jFeature}.fitResults = fR;
                
                % Update feature heirarchy to ensure all dependencies are acknowledged
                mainFeature.updateSubFeatures( );

                disp( sprintf('                            Exit Flag = %d', fR.exitflag) ); 
            end

        end
        % }}}
        % fitGlobal {{{
        function fitInfo = fitGlobal( obj, Image2Fit, parameters, addRem)

            if nargin == 6
                fitScope = ['globum_' addRem];
            else
                fitScope = 'global';
            end

            disp('                    Global Fitting : '); 
%             disp('                        Initiated'); 

            mainFeature = obj.featureList{ parameters.channelIdx, parameters.time};

            % Prepare for Fit (specialized)
            [ fitProblem, fitInfo] = prepareFit( obj, mainFeature, Image2Fit, parameters, fitScope);

            % Check if mainFeature should be fit. If not, then return out of method
            if isempty( mainFeature.featureList)
                fitInfo.featureMain = mainFeature;
                fprintf('                       Skipping fitting for %s\n', mainFeature.type);  
                return
            end

            % Run Global Fit
            [fR.vfit,fR.resnorm,fR.residual,fR.exitflag,~,~,fR.jacobian] = lsqnonlin( fitProblem ); 
            
            % Get the errors from the jacobian
            ci = nlparci( fR.vfit, fR.residual, 'jacobian', fR.jacobian, 'alpha', 0.07); fR.vfitError = abs( (ci(:,2)-ci(:,1) ) / 2 )';

            % Store fit Results
            fitInfo.fitResults = fR;

            disp( sprintf('                        Exit Flag = %d', fR.exitflag) ); 

        end
        % }}}
        % fitFeatureNumber {{{
        function fitInfo = fitFeatureNumber( obj, Image2Fit, parameters, fitInfo)

            % check if features should be added and removed
            if obj.settings.flags.skipOptimizeNumber
                return
            end

            cTime = parameters.time;
            cChannel = parameters.channelIdx;

            % check if mainFeature should be fit. If not, then return out of method
            if isempty( obj.featureList{ cChannel, cTime} )
                fprintf('                       Skipping fitting for %s\n', mainFeature.type);  
                return
            end

            % We will iteratively add and remove features to find the optimum number
            p = 1e-10;
%             alpha = 0.05;

            % Book-keeping
            fitInfo.fitScope = 'globum';
            fitInfo.Org = fitInfo;
            fitInfo.Old = fitInfo.Org;
            mainFeature.Org = obj.featureList{ cChannel, cTime};
            mainFeature.Old = copyDeep( mainFeature.Org);

            % Add Missing Features {{{
            continueAdd = 1;
            while continueAdd

                refImage = abs( Image2Fit - Cell.simulateCellWithVec( fitInfo.Old.fitResults.vfit, fitInfo.Old) );
               
                % Create a deep copy of the new  main feature
                mainFeature.New = mainFeature.Old.copyDeep(); 

                % Add Features
                mainFeature.New = mainFeature.New.addSubFeatures( refImage);
                
                % Run a global fit
                obj.featureList{ cChannel, cTime} = mainFeature.New;
                disp( sprintf('                - Add N = %d ---> %d ', mainFeature.Old.getSubFeatureNumber(), mainFeature.Old.getSubFeatureNumber()+1 ) )
                fitInfo.New = fitGlobal( obj, Image2Fit, parameters, 'add');
                
                im1_raw = mainFeature.Old.imageSim;
                im2_raw = mainFeature.New.imageSim;
                im1 = mainFeature.Old.simulateFeature( size( Image2Fit) );
                im2 = mainFeature.New.simulateFeature( size( Image2Fit) );
                numP1 = length( fitInfo.Old.fitVecs.vec);
                numP2 = length( fitInfo.New.fitVecs.vec);
                
                % Do a statistical F test to compare the two models
                p = Cell.compareModels(im1_raw, im2_raw, im1, im2, numP1, numP2, Image2Fit, 'f');
%                 disp(p)
%                 error('stop')

                if p > fitInfo.Org.alpha
                    disp(sprintf('                    N = %d (p = %.3f, alpha = %.3f) . Keep old model.', mainFeature.Old.getSubFeatureNumber(), p, fitInfo.Org.alpha ) )
%                     disp(sprintf('                    p = %.3f, alpha = %.3f, R1 = %.2f, R2 = %.2f', p, alpha, ss1, ss2 ) )
                    obj.featureList{cChannel,cTime} = mainFeature.Old;
                    continueAdd = 0;
                elseif p < fitInfo.Org.alpha
                    continueAdd = 1;
                    disp(sprintf('                    N = %d (p = %.3f, alpha = %.3f) . Keep improving model.', mainFeature.New.getSubFeatureNumber(), p, fitInfo.Org.alpha) )
%                     disp(sprintf('                    p = %.3f, alpha = %.3f, R1 = %.2f, R2 = %.2f', p, alpha, ss1, ss2 ) )                    
                    fitInfo.Old = fitInfo.New;
                    obj.featureList{cChannel,cTime} = mainFeature.New;
                    mainFeature.Old = mainFeature.New;
                end

            end
            % }}}

            % Remove Redundant Features {{{
            continueRemove = 1;
            while continueRemove
               
                refImage = abs( Image2Fit - Cell.simulateCellWithVec( fitInfo.Old.fitResults.vfit, fitInfo.Old) );

                % Create a deep copy of the new  main feature
                mainFeature.New = copyDeep( mainFeature.Old); 
                
                % Remove Features
                [ mainFeature.New, successRemove ] = mainFeature.New.removeSubFeatures( refImage);
                
                % A removal could fail if there are no features left to remove
                if successRemove

                    % Run a global fit
                    disp( sprintf('                - Rem N = %d ---> %d ', mainFeature.Old.getSubFeatureNumber(), mainFeature.Old.getSubFeatureNumber()-1 ) )
                    obj.featureList{ cChannel, cTime} = mainFeature.New;
                    if mainFeature.New.getSubFeatureNumber() > 0
                        fitInfo.New = fitGlobal( obj, Image2Fit, parameters, 'remove');
                    end
                    
                    im1_raw = mainFeature.Old.imageSim;
                    im2_raw = mainFeature.New.imageSim;
                    im1 = mainFeature.Old.simulateFeature( size( Image2Fit) );
                    im2 = mainFeature.New.simulateFeature( size( Image2Fit) );
                    numP1 = length( fitInfo.Old.fitVecs.vec);
                    numP2 = length( fitInfo.New.fitVecs.vec);
                    
                    % Do a statistical F test to compare the two models
                    p = Cell.compareModels(im2_raw, im1_raw, im2, im1, numP2, numP1, Image2Fit, 'f');
                    
                    if p < fitInfo.Org.alpha
                        disp(sprintf('                    N = %d (p = %.3f, alpha = %.3f) . Keep old model.', mainFeature.Old.getSubFeatureNumber(), p, fitInfo.Org.alpha) )
%                         disp(sprintf('                    p = %.3f, alpha = %.3f, R1 = %.2f, R2 = %.2f', p, alpha, ss1, ss2 ) )                    
                        obj.featureList{cChannel,cTime} = mainFeature.Old;
                        continueRemove = 0;
                    elseif p > fitInfo.Org.alpha
                        continueRemove = 1;
                        disp(sprintf('                    N = %d (p = %.3f, alpha = %.3f) . Keep improving model.', mainFeature.New.getSubFeatureNumber(), p, fitInfo.Org.alpha ) )
%                         disp(sprintf('                    p = %.3f, alpha = %.3f, R1 = %.2f, R2 = %.2f', p, alpha, ss1, ss2 ) )                    
                        fitInfo.Old = fitInfo.New;
                        obj.featureList{cChannel,cTime} = mainFeature.New;
                        mainFeature.Old = mainFeature.New;
                    end

                else
                    continueRemove = 0;
                    disp( sprintf('                - Rem N = %d ---> %d ', mainFeature.Old.getSubFeatureNumber(), mainFeature.Old.getSubFeatureNumber()-1 ) )
                    disp( sprintf('                    No more removable features', mainFeature.Old.getSubFeatureNumber() ) )
                    disp( sprintf('                    N = %d , Keep old model.', mainFeature.Old.getSubFeatureNumber() ) )
                end

            end
            % }}}

            disp( sprintf('                - FINAL N =  %d', mainFeature.Old.getSubFeatureNumber() ) )
            fitInfo = fitInfo.Old;
            obj.featureList{ cChannel, cTime} = mainFeature.Old;

        end
        % }}}
        % prepareFit {{{
        function [ fitProblem, fitInfo ] = prepareFit( obj, featureMain, imageOrg, parameters, fitScope, cFeature)
            
            % Input Checks
            if length( featureMain) ~= 1
                error('prep_fitLocal : there should only be a single main feature'), end
            
            if ~strcmp( fitScope, 'local') && ~strcmp( fitScope, 'global') && ~strcmp( fitScope, 'globum_add') && ~strcmp( fitScope, 'globum_remove')
                error('prepareFit : input argument fitScope must be either ''local'' or ''global'' '); end

            if strcmp( fitScope, 'local') && nargin < 6
                error('prepareFit : for a local fit, feature number argument must be provided'); end

            % The fit will be 3-dimensional
            imageOrg = im2double( imageOrg); 
        
            % Set optimization options based on configuration flags
%             opts = Cell.setOptimOpts( obj.settings.flags);

            % If no features are present to be fitted, assign the essentials and return out
            if isempty( featureMain.featureList)
                fitInfo.channel = parameters.channelTrue;
                fitInfo.time = parameters.time;
                fitInfo.fitScope = fitScope;
                fitInfo.saveDirectory = parameters.saveDirectory;
                fitProblem = [];
                return
            end

            % Obtain the feature to fit, its fit vector and its labels for lsqnonlin and parsing
            if strcmp(fitScope, 'local')
                [ fitVec, fitLabels, fitObj] = getVecLocal( featureMain);
                fitVec = fitVec{ cFeature};
                fitLabels = fitLabels{ cFeature};
                fitObj = fitObj{ cFeature};
            elseif strcmp(fitScope, 'global') || strcmp(fitScope, 'globum_add') || strcmp(fitScope, 'globum_remove')
                [ fitVec, fitLabels ] = getVec( featureMain); 
                fitObj = featureMain; cFeature = 1; 
                fitInfo.Nnew = featureMain.getSubFeatureNumber();
            end

            % Record old number of features for optimization of feature number
            if strcmp(fitScope, 'globum_add')
                fitInfo.Nold = fitInfo.Nnew - 1;
            elseif strcmp(fitScope, 'globum_remove')
                fitInfo.Nold = fitInfo.Nnew + 1;
            end
            % Get upper and lower bounds of parameters for passing on to lsqnonlin (restricts the parameter space)
            fitVecs = Cell.getUpperLowerBounds( fitVec, fitLabels, imageOrg);

            % Scale the parameters to vary the speed of exploration in the parameter space
            if obj.settings.flags.fitExploreSpeed
                speedVec = Cell.getExplorationSpeedVector( fitLabels);
                fitVecs.vec = fitVecs.vec ./ speedVec;
                fitVecs.ub = fitVecs.ub ./ speedVec;
                fitVecs.lb = fitVecs.lb ./ speedVec;
                fitInfo.speedVec = speedVec;
            end

            % Set up parameters for fit
            fitInfo.featureMain = featureMain;
            fitInfo.featureCurrent = fitObj;
            fitInfo.featureIndex = cFeature;
            fitInfo.fitVecs = fitVecs;
            fitInfo.mask = logical( imageOrg);
            fitInfo.image = imageOrg;
            fitInfo.numVoxels = size( imageOrg);
            fitInfo.channel = parameters.channelTrue;
            fitInfo.time = parameters.time;
            fitInfo.fitScope = fitScope;
            fitInfo.saveDirectory = parameters.saveDirectory;
            fitInfo.alpha = 0.05;
            % Find voxel indices for gaussian feature computation

            % make the error function for lsqnonlin
            f = Cell.makeErrorFcn( imageOrg, fitInfo );

            % Crate Optimization Problem

            % Set Optim Options {{{
            flags = obj.settings.flags;
            if flags.debug 

%                                     'UseParallel', true, ...
                opts = optimoptions( @lsqnonlin, ...
                                    'MaxFunEvals', 2000, ...
                                    'OptimalityTolerance', 1e-12, ...
                                    'MaxIter', 20, ...
                                    'TolFun', 1e-9, ...
                                    'FiniteDifferenceStepSize', 1e-2, ...
                                    'FiniteDifferenceType', 'central', ...
                                    'StepTolerance', 1e-5, ...
                                    'display', 'iter', ... 
                                    'OutputFcn', @plotFit );

            elseif flags.graphRealTime

                opts = optimoptions( @lsqnonlin, ...
                                    'MaxFunEvals', 2000, ...
                                    'OptimalityTolerance', 1e-12, ...
                                    'MaxIter', 10, ...
                                    'TolFun', 1e-9, ...
                                    'FiniteDifferenceStepSize', 1e-2, ...
                                    'FiniteDifferenceType', 'central', ...
                                    'StepTolerance', 1e-5, ...
                                    'display', 'off', ... 
                                    'OutputFcn', @plotFit );

            else 
                opts = optimoptions( @lsqnonlin, ...
                                    'MaxFunEvals', 2000, ...
                                    'OptimalityTolerance', 1e-12, ...
                                    'MaxIter', 30, ...
                                    'TolFun', 1e-9, ...
                                    'FiniteDifferenceStepSize', 1e-2, ...
                                    'FiniteDifferenceType', 'central', ...
                                    'StepTolerance', 1e-5, ...
                                    'display', 'iter', ...
                                    'UseParallel', false);
            end
            % }}}

            fitProblem = createOptimProblem( 'lsqnonlin', 'objective', f, 'x0', fitVecs.vec, 'ub', fitVecs.ub, 'lb', fitVecs.lb, 'options', opts);

            % Create handle for plotting function
            function stop = plotFit( x, optimV, state)
                stop = Cell.plot_midFit(x, optimV, state, fitInfo);
            end

        end
        % }}}
        
        % }}}
        
        % playMovie {{{
        function playMovie( obj, cChannels, cColors)
            % Plays the image movie of the Cell.
            % cChannels can either be a scalar or a vector of upto 3 elements specifying the channels to play
            % cColor can either be a string or a string array of upto 3 elements specifying the color of each of the channels. the size of cColors must match the size of cChannels

            if length( cChannels) > 3
                error('playMovie: cChannels cannot contain more than 3 elements')
            end

            if nargin == 3 && length( cChannels) ~= length( cColors)
                error('playMovie: length of cChannels must match length of cColors')
            end
            
            img = max( obj.image, [], 3);
            img = squeeze( img);
            img = permute( img, [ 1 2 4 3] ); % XYCT

            imStack = zeros( size(img, 1), size(img,2), 3, size(img, 4) );
            rgb = {'R', 'G', 'B'};

            if nargin < 3
                if length( cChannels) == 1
                    % make greyscale image
                    imStack(:,:,1,:) = mat2gray( img(:,:,cChannels, :) );
                    imStack(:,:,2,:) = mat2gray( img(:,:,cChannels, :) );
                    imStack(:,:,3,:) = mat2gray( img(:,:,cChannels, :) );
                else
                    for jC = 1 : length(cChannels)
                        imStack(:,:,jC,:) = mat2gray( img(:,:,cChannels(jC), :) );
                    end
                end

            elseif nargin == 3
            
                % Find cChannel for red color
                cR = cChannels( find( strcmp(cColors, 'R') ) );
                if ~isempty( cR)
                    imStack(:,:,1,:) = mat2gray( img(:,:,cR, :) );
                end

                % Find cChannel for green color
                cG = cChannels( find( strcmp(cColors, 'G') ) );
                if ~isempty( cG)
                    imStack(:,:,2,:) = mat2gray( img(:,:,cG, :) );
                end
                
                % Find cChannel for blue color
                cB = cChannels( find( strcmp(cColors, 'B') ) );
                if ~isempty( cB)
                    imStack(:,:,3,:) = mat2gray( img(:,:,cB, :) );
                end
            end

            fps = 10;
            implay( imStack, fps);

        end
        % }}}
        
        % simulateCell {{{
        function imageOut = simulateCell( obj, cChannel, cTime)
            % Simulates XYZ image of a cell with a mainFeature
            
            mainFeature = obj.featureList{cChannel, cTime};

            imageOut = 0*mainFeature.image;

            % simulate the main feature 
            imageOut = mainFeature.simulateFeature( imageOut);

        end
        % }}}
        
        % updateFeatureMap {{{
        function obj = updateFeatureMap( obj, channel, time)
           
            % Initialize the counter. Get all the keys from featureIdxList
            % and add 1
            mapGlobal = obj.featureMap;
            counter = max( cell2mat( keys( mapGlobal) ) );
            if isempty( counter)
                counter = 0;
            end
            
            % Get local book, its keys, and its values
            mapLocal = obj.featureList{channel,time}.featureMap;
            keysLocal = keys( mapLocal);
            valsLocal = values( mapLocal);
            
            % update local keys by counter
            % add values to the global map
            if ~isempty( keysLocal);
                keysLocalNew = cell2mat(keysLocal) + counter;
            end
            remove( mapLocal, keysLocal);
            for jKey = cell2mat(keysLocal)
                % add to local map
                mapLocal( jKey) = valsLocal{jKey};
                % add to global map
                 % Now with the updated local hashmap, we will prepend the
                % channel and time to the values and add the new keys and new
                % values to the global hashmap
                mapGlobal( jKey) = [ channel time valsLocal{jKey}];
            end
            
            % ----------------------------------------------------
            % FORCE SUBFEATURE ID UPDATES WHEN MAPLOCAL KEYS CHANGE
            % ----------------------------------------------------
            
        end
        % }}} 

    end
    
    methods ( Static = true , Access = public )
        
       % Common Functions {{{
       
        % img2double {{{
        function imgOut = img2double( imgIn)

            mx = double( max( imgIn(:) ) );
            mn = double( min( imgIn(:) ) );
            dif = mx-mn;
            imgOut = mn + dif*im2double( imgIn);

        end
        % }}}
        % img2uint16 {{{
        function imgOut = img2uint16( imgIn)

            mx = double( max( imgIn(:) ) );
            mn = double( min( imgIn(:) ) );
            dif = mx-mn;
            imgOut = mn + dif*im2double( imgIn);

        end
        % }}}
        % findAmplitudeAlongLine {{{
        function amplitude = findAmplitudeAlongLine( imageIn, startPos, endPos)

            dim = length( size( imageIn) );
            if dim ~= length( startPos) || dim ~= length( endPos)
                error( 'findAmplitudeAlongLine: dimensions of input arguments does not match') 
            end

            imageIn = im2double( imageIn);

            LineLength = norm( abs(endPos - startPos) );

            numPoints = round( LineLength*5);
            if dim == 2
                xq = linspace( startPos(1), endPos(1), numPoints);
                yq = linspace( startPos(2), endPos(2), numPoints);
                amplitude = interp2( imageIn, yq, xq, 'linear');

            elseif dim == 3
                xq = linspace( startPos(1), endPos(1), numPoints);
                yq = linspace( startPos(2), endPos(2), numPoints);
                zq = linspace( startPos(3), endPos(3), numPoints);
                amplitude = interp3( imageIn, xq, yq, zq, 'linear');

            else
                error('findAmplitudeAlongLine: dimensionality of data must be 2D or 3D')
            end

        end
        % }}}
        % radIntegrate2D {{{
        function [phiIntensity, phiValues] = radIntegrate2D( imageIn, startPoint)

            if length( size(imageIn) ) ~= 2
                error( 'radIntegrate2D : input image must be 2-dimensional')
            end

            % Radially integrate around startPoint 
            PhiStep = deg2rad(2);
            rmin = 1;
            rmax = 10;
            numVoxelsX = size( imageIn, 2);
            numVoxelsY = size( imageIn, 1);
            x0 = startPoint(1);
            y0 = startPoint(2);

            % Define the step sizes to use for the angular sweep in phi. 
            PhiStep = deg2rad(2);
            rmin = 1;
            rmax = 10;
            phiValues = 0 : PhiStep : 2*pi - PhiStep;

            % Pre-allocate array for storing intensity data with varying phi.
            phiIntensity = zeros( 1, length( phiValues) );

            method = 'interp2';
            switch method

            % method: interp2 {{{
            case 'interp2'
                for jPhi = 1 : length( phiValues)

                    % Find the coords along this line
                    X1 = x0 + ( [rmin:0.5:rmax] .* cos( phiValues( jPhi) ) );
                    Y1 = y0 + ( [rmin:0.5:rmax] .* sin( phiValues( jPhi) ) ); 

                    % Identify coords that exceed the size of the image, and remove them
                    rmIdx = union( find( X1 < 1 | X1 > numVoxelsX) , find( Y1 < 1 | Y1 > numVoxelsY) );
                    X1( rmIdx) = [];
                    Y1( rmIdx) = [];
                    
                    % line Intensity along the line
                    lineInt = interp2( im2double( imageIn), X1, Y1, 'linear');
                    phiIntensity( jPhi) = mean( lineInt);
                    
                end
            % }}}

            % method: convexhull {{{
            case 'convexhull'
                for jPhi = 1 : length( phiValues)
                    % Extract the relevant angles for drawing a non-zero angle
                    phi1 = phiValues( jPhi) - PhiStep;
                    phi2 = phiValues( jPhi) + PhiStep;

                    % Find the coordinates of the extreme points on the arc length
                    % of this pizza slice.
                    X1 = x0 + ( rmax * cos( [ phi1, phi2]) );
                    Y1 = y0 + ( rmax * sin( [ phi1, phi2]) ); % -ve because matlab orientation is reversed( origin is at top left)

                    % find minimum points at distance rmin away
                    X2 = x0 + ( rmin * cos( [ phi1, phi2]) );
                    Y2 = y0 + ( rmin * sin( [ phi1, phi2]) ); 

                    X = [ X1, X2];
                    Y = [ Y1, Y2];

                    xRd = []; yRd = [];
                    for jPt = 1 : length(X)
                        [~, xidx] = min( abs( xVecNew - round( X(jPt), 1) ) );
                        [~, yidx] = min( abs( xVecNew - round( Y(jPt), 1) ) );
                        xRd = [ xRd, xidx ];
                        yRd = [ yRd, yidx ];
                    end

                    % Initialize the mask. We'll turn on the (X,Y) pixels and draw a convex
                    % hull around it.
                    imMask = zeros( size( imSub) );
                    imMask( sub2ind( size( imSub), yRd, xRd ) ) = 1;

                    % Draw a convex hull and uncover the pixel values of interest
                    imMask = bwconvhull( imMask);

                    imMasked = imMask .* imSub;
            %         imshowpair(imSub, imMask); pause(0.25)
                    
                    % Sum up values and store them
                    phiIntensity( jPhi) = sum( imMasked(:) ) / sum(imMask(:));
                    
                end
                % }}}

            end

            % We smooth the angular intensity to enable easy peak finding
            phiIntensity = imgaussfilt( phiIntensity, 3);
           
            % find idx of minimum intensity and shift the array so no peaks are on the edge
            [~, idxShift] = min( phiIntensity);
            phiIntensity = circshift( phiIntensity, -idxShift);
            phiValues( 1:idxShift) = phiValues(1:idxShift) + 2*pi;
            phiValues = circshift( phiValues, -idxShift);
%             figure; plot( phiValues, phiIntensity, 'r-')

        end
        % }}}
        % drawGaussianPoint3D {{{
        function imagePoint = drawGaussianPoint3D( pos, sigma, imageIn, idx, x, y, z)
            % Draws a gaussian point using an analytical framework

            % Check that everything is in 3D
            if numel( size(imageIn) ) ~= 3 || length(pos)~=3 || length(sigma)~=3
                error('drawGaussianLine3D: all data must be 3-dimensional')
            end

            if nargin < 4
                idx = 1 : numel( imageIn);
            end
            if nargin < 7
                [y, x, z] = ind2sub( size( imageIn), idx);
            end

            x0 = pos(1); y0 = pos(2); z0 = pos(3);
            sx = sigma(1); sy = sigma(2); sz = sigma(3);

            % Amplitudes:
            ExpX = exp( -0.5*( (x-x0)./sx).^2 );
            ExpY = exp( -0.5*( (y-y0)./sy).^2 );
            ExpZ = exp( -0.5*( (z-z0)./sz).^2 );
            IntValues = ExpX .* ExpY .* ExpZ;
           
            functionCheck = 0;
            if functionCheck
                disp( sprintf('Total Voxels = %d', numel( IntValues) ) )
                disp( sprintf('Nan Voxels = %d', sum( isnan( IntValues(:) ) ) ) )
                disp( sprintf('Inf Voxels = %d', sum( isinf( IntValues(:) ) ) ) )
            end
            IntValues( isnan( IntValues) ) = min( IntValues(:) );
            IntValues( IntValues == Inf) = min( IntValues(:) );
            
            % Initiliaze the volume and set the appropriate indices to these values
            imagePoint = 0*imageIn;
            imagePoint( idx) = IntValues;

            if functionCheck
                dispImg( imagePoint);
            end
                    
        end
        % }}}
        % drawGaussianLine3D {{{
        function imageLine = drawGaussianLine3D( startPos, endPos, sigma, imageIn, idx, x, y, z)
            % Draws a gaussian straight line using an analytical framework

            % Check that everything is in 3D
            if numel( size(imageIn) ) ~= 3 || length(startPos)~=3 || length(endPos)~=3 || length(sigma)~=3
                error('drawGaussianLine3D: all data must be 3-dimensional')
            end

            if nargin < 5
                idx = 1 : numel( imageIn);
            end
            if nargin < 8
                [y, x, z] = ind2sub( size( imageIn), idx);
            end

            x0 = startPos(1); y0 = startPos(2); z0 = startPos(3);
            x1 = endPos(1); y1 = endPos(2); z1 = endPos(3);
            sx = sigma(1); sy = sigma(2); sz = sigma(3);

            % First lets parameterize this line with a parameter t. We'll find the
            % slopes of the line in each spatial dimension along with any offset
            % At t = 0, x = x0, y = y0, z = z0
            % At t = 1, x = x1, y = y1, z = z1
            t0 = 0;
            t1 = 1;
            
            % find slopes
            mx = (x1 - x0) / ( t1 - t0);
            my = (y1 - y0) / ( t1 - t0);
            mz = (z1 - z0) / ( t1 - t0);
            
            % find offsets, these are just the start points at t = 0;
            cx = x0;
            cy = y0;
            cz = z0;
            
            % We will integrate a 3D gaussian with respect to the parameter t going
            % from 0 to 1.
            
            % We've done the analytical integration in mathematica and will use the
            % result here:
            
            % Amplitudes:
            AmpExp = - ( 1 / sqrt(mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2 ) );
            
            AmpErf = sqrt( pi/2 ) * sx * sy * sz;
            
            % Exponential factors:
            ExpDenom = 2 * ( mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2);
            
            Exp1 = ( - ( cx^2 * mz^2 * sy^2 + cz^2 * (my^2 * sx^2 + mx^2 * sy^2) + ...
                cx^2 * my^2 * sz^2 + cy^2 * (mz^2 * sx^2 + mx^2 * sz^2) ) / ExpDenom );
            
            Exp2 = ( - ( - 2 * cx * mz^2 *sy^2 * x - 2 * cx * my^2 * sz^2 * x + ...
                mz^2 * sy^2 * x.^2 + my^2 * sz^2 * x.^2 + 2 * cx * mx * my * sz^2 * y - ...
                2 * mx * my * sz^2 * x .* y + mz^2 * sx^2 * y.^2 + mx^2 * sz^2 * y.^2 ) / ExpDenom  );
            
            Exp3 = ( - ( 2 * cx * mx * mz * sy^2 * z - 2 * mx * mz * sy^2 * x .* z - ...
                2 * my * mz * sx^2 * y .* z + my^2 * sx^2 * z.^2 + mx^2 * sy^2 * z.^2  ) / ExpDenom );
            
            Exp4 = ( - ( - 2 * cz * ( cy * my * mz * sx^2 + cx * mx * mz * sy^2 - ...
                mx * mz * sy^2 * x - my * mz * sx^2 * y + my^2 * sx^2 * z + mx^2 * sy^2 * z )  ) / ExpDenom );
            
            Exp5 = ( - ( - 2 * cy * (cx *mx *my *sz^2 - mx * my * sz^2 * x + ...
                mx^2 * sz^2 * y + mz * sx^2 * (mz * y - my * z)) ) / ExpDenom );
           
            SumExp = Exp1 + Exp2 + Exp3 + Exp4 + Exp5;
            
        %     Exp5( Exp5 == Inf) = realmax;
            
            % Erf factors;
            ErfDenom = ( sqrt(2) * sx * sy * sz * sqrt( mz^2 * sx^2 * sy^2 + (my^2 * sx^2 + mx^2 * sy^2) * sz^2 ) );
            
            Erf1 = erf( ( cz * mz * sx^2 * sy^2 + cy * my * sx^2 * sz^2 + cx * mx * sy^2 * sz^2 - ...
                mx * sy^2 * sz^2 * x - my * sx^2 * sz^2 * y - mz * sx^2 * sy^2 * z ) / ErfDenom );
            
            Erf2 = erf( ( cz * mz * sx^2 * sy^2 + mz^2 * sx^2 * sy^2 + ...
                sz^2 * (cy * my * sx^2 + my^2 * sx^2 + mx * sy^2 * (cx + mx - x) - my * sx^2 * y) - ...
                mz * sx^2 * sy^2 * z ) / ErfDenom );

            % Find the intensity values for these provided query x,y,z coordinates
            IntValues = AmpExp .*exp( SumExp) .* ...
                AmpErf .* ( Erf1 - Erf2 );
           
            functionCheck = 0;
            if functionCheck
                disp( sprintf('Total Voxels = %d', numel( IntValues) ) )
                disp( sprintf('Nan Voxels = %d', sum( isnan( IntValues(:) ) ) ) )
                disp( sprintf('Inf Voxels = %d', sum( isinf( IntValues(:) ) ) ) )
                figure; imagesc( max( imageLine, [], 3) );
                figure; imagesc( max( imageLine, [], 3) );
            end

            IntValues( isnan( IntValues) ) = min( IntValues(:) );
            IntValues( IntValues == Inf) = min( IntValues(:) );
            
            % Initiliaze the volume and set the appropriate indices to these values
            imageLine = 0 * imageIn;
            imageLine( idx) = IntValues;

        end
        % }}}
        % simulateCellWithVec {{{
        function imageOut = simulateCellWithVec( vec, fitInfo)
            % Simulates the Features in the cell
            
            % Apply the speedExploration vector to transform to real parameters
            if isfield( fitInfo, 'speedVec')
                vec = vec .* fitInfo.speedVec;
            end

            % Get current feature
            featureCurrent = fitInfo.featureCurrent;
            featureCurrent.absorbVec( vec, fitInfo.fitVecs.labels );
            
            % Simulate the feature
            imageOut = fitInfo.featureMain.simulateAll( 0*fitInfo.image, featureCurrent.ID);

        end
        % }}}
            % checkFrameViability {{{
            function status = checkFrameViability( Image2Fit, jTime);

                status = 0;
                % Check if frame is all super-bright( i.e. voxel variance is close to zero)
                maxVox = max( Image2Fit(:) );
                minVox = min( Image2Fit(:) );

                if maxVox == minVox
                    status = 1
                    disp( sprintf( 'Unviable frame %d, Skipping it...', jTime) )
                end

            end
            % }}}
            
       % }}}
       
        % Feature Estimation {{{
        % find_manyLines {{{
        function [Lines, ax] = find_manyLines( imageIn, startPoint, ax)

            imageIn = im2double( imageIn);
            image2D = max( imageIn, [], 3);

            if nargin < 3 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( imageIn, [], 3) ); colormap gray; axis equal; hold on;
            end

            % radially integrate in phi
            [ phiIntensity, phiValues] = Cell.radIntegrate2D( image2D, startPoint);
%             figure; plot( phiValues, phiIntensity, 'r-');

            % find peaks in Phi Intensity
            imVals = image2D( image2D ~= 0);
            mtBkg = median( imVals );
            minPkHeight = mtBkg + std( imVals );
            
            warning('off', 'signal:findpeaks:largeMinPeakHeight' )
            [ peakIntensity, peakPhi ] = findpeaks( phiIntensity, phiValues, 'MinPeakHeight', minPkHeight );
            peakPhi = mod( peakPhi, 2*pi);
            warning('on', 'signal:findpeaks:largeMinPeakHeight' )

            Lines = {};
            % Find length of lines
            for jLine = 1 : length( peakPhi)

                if nargout > 1
                    [ cLine, ax] = Cell.find_singleLine3D( imageIn, startPoint, peakPhi( jLine), ax );
                else
                    cLine = Cell.find_singleLine3D( imageIn, startPoint, peakPhi( jLine) );
                end

                if ~isempty( cLine)
                    Lines = { Lines{:}, cLine}; 
                end

            end
%             error('stop here')
        end
        % }}}
        % find_singleLine3D {{{
        function [Line, ax] = find_singleLine3D( imageIn, startPoint, startAngle, ax)
            % This works with both 3D and a 2D input input. The outputs vary for 2D and 3d inputs

            if nargin < 4 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( image2D, [], 3) ); colormap gray; axis equal; hold on;
            end

            dim = length( size( imageIn) );
            [image2D, idx3Max] = max( imageIn, [], 3);
            
            startPoint2D = startPoint(1:2);

            if nargout > 1
                [Line, ax] = Cell.find_singleLine2D( image2D, startPoint2D, startAngle, ax);
            else
                Line = Cell.find_singleLine2D( image2D, startPoint2D, startAngle);
            end

            if isempty( Line)
                return
            end

            if dim == 3
                % Find Z position for this Line, find theta angle for this Line

                % Approach : from the computed MIP, we also have access to the respective
                % z-indices that correspond to each maximum value. For each guessed
                % microtubule, we'll point it in its direction and get the local indices to
                % the nearest pixel. With these indices we can compute what theta must be.

                % define start radius for theta determination
                rmin = 1;

                % Define the radia up to which you'll check the z-positions
                RadVec = rmin : Line.length;
                ThetaVec = 0*RadVec;

                % for every possible point on the Line before it ends
                % final coords
                x1 = startPoint(1) + RadVec * cos( startAngle );
                y1 = startPoint(2) + RadVec * sin( startAngle );

                % get the z coord of the maxidx at this position, and calculate theta based on this value.
                z0 = startPoint(3);
%                 z0 = idx3Max( startPoint(2), startPoint(1) );
                z1 = interp2( idx3Max, x1, y1, 'Linear');
                zdist = abs( z1 - z0);
                theta = atan( abs( z1 - z0) ./ RadVec);
               
                % Based on how theta is defined, starting at the positive z-axis
                % and sweeping pi radians until it hits the negative z-axis, we'll
                % find theta depending on the location of z1 above or below the SPB
                % in the z-dimension
                theta( z1 >= z0) = pi/2 - theta( z1 >= z0);
                theta( z1 < z0) = pi/2 + theta( z1 < z0);
                Line.theta = mean(theta);
                
                z1 = z0 + Line.length*cos( Line.theta);
                Line.startPosition = [ Line.startPosition , z0 ];
                Line.endPosition = [ Line.endPosition , z1 ];
            end

        end
        % }}}
        % find_singleLine2D {{{
        function [Line, ax] = find_singleLine2D( image2D, startPoint, startAngle, ax)

            % Input Checks
            if length( startPoint) ~= 2
                error('find_singleLine2D: input argument startPointmust be of length 2')
            end

            if nargin < 4 && nargout > 1
                h = figure('Visible', 'on');
                ax = axes;
                imagesc( max( image2D, [], 3) ); colormap gray; axis equal; hold on;
            end

            Lmin = 5;

            % This will look at radial amplitude alon the Line to find when the intensity falls below a threshold
            threshIntensity = median( image2D( image2D ~= 0) ) + 1.5*std( image2D( image2D ~= 0) );
%             median( image2D( image2D ~= 0) ) 
%             std( image2D( image2D ~= 0) )
%             error('stop here')

            % Define what radia to use:
            rmin = 1;
            rmax = 20;
            radialValues = rmin : 0.5 : rmax;
            radialIntensity = 0* radialValues;
            x0 = startPoint( 1);
            y0 = startPoint( 2);
                
            for jRad = 1 : length( radialValues)

                jLen = radialValues( jRad);
                % Extract the relevant angles for drawing a non-zero angle

                % Find the coordinates of the extreme points on the arc length of this.
                X = x0 + ( jLen * cos( startAngle) );
                Y = y0 + ( jLen * sin( startAngle) );

                % Get the intensity value at this position
                radialIntensity( jRad) = interp2( image2D, X, Y, 'Linear');
                
            end

            % get the length value of this Line by finding the length along
            % the Line when the intensity drops below a certain threshold.
%             radialIntensity = smooth( radialIntensity, 3);

%             figure; plot(radialValues, radialIntensity); 
        
            [ ~, IdMax] = find( radialIntensity < threshIntensity, 1, 'first');
            if all(radialIntensity > threshIntensity), IdMax = length(radialIntensity); end

            if isempty( IdMax) || radialValues( IdMax) < Lmin, 
                Line = [];
                return
            else 
                % Store the corresponding properties for the Line 
                Line.length = radialValues( IdMax);
                Line.phi = startAngle;
                Line.startPosition = startPoint;
                Line.endPosition = startPoint + Line.length * [cos(startAngle), sin(startAngle)];
                if nargout > 1
                    line( [Line.startPosition(1), Line.endPosition(1)], [Line.startPosition(2), Line.endPosition(2)], 'Color', 'r', 'LineWidth', 3); 
                end

            end

        end
        % }}}
        % }}}
        
        % Feature Fitting {{{

        % getUpperLowerBounds {{{
        function fitVecs = getUpperLowerBounds( vec, vecLabels, image)
            % Get upper and lower bounds for fitting
            
            ub = vec; lb = vec;
            minVox = 1;
            maxVox = size( image);
            estBkg = median( image( image> 0) );
            maxBkg = max( image(:) );
            minBkg = min( image(:) );
            maxAmp = max( image(:) );
            minAmp = 0; % min SnR is half of max SnR of image
            minSig = [1.2, 1.2, 1.0];
            maxSig = [2.0, 2.0, 2.0];
            

            % find the correct label in vecLabels, and start to place in the correct bounds in the correct places


            % Find index of parameters 
            idxAmp = find( ~cellfun( @isempty, strfind( vecLabels, 'amplitude') ) );
            idxSig = find( ~cellfun( @isempty, strfind( vecLabels, 'sigma') ) );
            idxP0 = find( ~cellfun( @isempty, strfind( vecLabels, 'startPosition') ) );
            idxP1 = find( ~cellfun( @isempty, strfind( vecLabels, 'endPosition') ) );
            idxP = find( ~cellfun( @isempty, strfind( vecLabels, 'position') ) );
            idxBkg = find( ~cellfun( @isempty, strfind( vecLabels, 'background') ) );

            % determine amplitude bounds
            % halfway between background and max value?
            %minAmp = 0.5*(maxBkg - vec(1) );
            %idxUp = vec(idxAmp) < minAmp;
            %vec( idxAmp(idxUp) ) = minAmp;

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
                ub( idxP ) = maxBkg;
                lb( idxP ) = minBkg; end

            if any( lb > ub) || any(ub < lb) || any( vec < lb) || any(vec > ub)
                badIdx = unique([ find( lb > ub),  find(ub < lb) , find( vec < lb) , find(vec > ub)]);
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
        % makeErrorFcn {{{
        function errVal = makeErrorFcn( ImageFit, fitInfo)

            errVal = @ErrorFcn;

            function err = ErrorFcn( p)
                    
                err = Cell.simulateCellWithVec( p, fitInfo) - ImageFit;
                err = err(:);

            end

        end
        % }}}
        % compareModels {{{
        function vargout = compareModels(im1_raw, im2_raw, im1, im2, numP1, numP2, imOrig, test)
            
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
        
        % }}}
            
        % Graphics {{{
        
        % plot_midFit {{{
        function stop = plot_midFit( vec, optim, state, fitInfo )
            % Plots during Lsqnonlin fitting
            % Inputs :  
            %   OPTIMVALUES: properties are funccount, iteration
            %   STATE: string representing current state of optimization engine. Values are init, iter, done 
            %   STOP: A boolean to stop the algorithm.
            %   All other inputs available to error function are also available. Currently: ImageOrg( original image), fitInfo and parameters

            stop = false;

            ImageOrg = fitInfo.image;
            nX = fitInfo.numVoxels(2); nY = fitInfo.numVoxels(1); nZ = fitInfo.numVoxels(3); 
            dim = length( size(ImageOrg) );
            image2D = max( ImageOrg, [], 3); intMin = min( ImageOrg(:) ); intMax = max( ImageOrg(:) );
            
            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

            time = fitInfo.time;
            channel = fitInfo.channel;
            fitScope = fitInfo.fitScope;
            if strcmp( fitScope, 'local')
                fName = sprintf('C%d_T%d_%s_feat%d_iter%d', channel, time, fitScope, fitInfo.featureIndex, optim.iteration );
                sName = sprintf('%s_feat%d_iter%d',  fitScope, fitInfo.featureIndex, optim.iteration );
            elseif strcmp( fitScope, 'global')
                fName = sprintf('C%d_T%d_%s_iter%d', channel, time, fitScope, optim.iteration );
                sName = sprintf('%s_iter%d', fitScope, optim.iteration );
            elseif strcmp( fitScope, 'globum_add')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            elseif strcmp( fitScope, 'globum_remove')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            end
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            graphicsVerbose = 1;
            if graphicsVerbose == 1
                %  Interactive Plotting {{{
                switch state
                    case 'init'

                        % Make Figure
                        figTitle = [];
                        figName = [];
                        figPath = [];
                        f = figure( figProps{:} );
                        set( f, 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope) );
                        ax = tight_subplot(2,3, 0.05);
                        delete( ax(6) );
                        drawnow
                        pause(1)                        

                        % Original Image
                        axes( ax(1) ); %                     subplot(231)
                        img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');
                        set( ax(1), 'Tag', 'ax1');

                        % Simulated Image
    %                     subplot(232)
                        axes( ax(2) ); 
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = imagesc( max(imageSim, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated'); set( img, 'Tag', 'imgSim');
                        set( ax(2), 'Tag', 'ax2');

                        % Residual Image
    %                     subplot(233)
                        axes( ax(3) ); 
                        img = imagesc( max( abs( ImageOrg - imageSim), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Residual'); set( img, 'Tag', 'imgRes');
                        set( ax(3), 'Tag', 'ax3');

                        % Graphical Features
    %                     subplot(234)
                        axes( ax(4) ); 
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(4), 'Tag', 'ax4');

                        % Residual vs Iteration
    %                     subplot(235)
                        axes( ax(5) ); 
                        xtickformat( '%.2f'); ytickformat( '%d')
                        plotRes = plot( optim.iteration, optim.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set( plotRes,'Tag','plotRes');
                        xlabel('Iteration','interp','none'); ylabel('Resnorm','interp','none');
                        title( sprintf('Best Resnorm : %g', optim.resnorm ),'interp','none');
                        set( gca, 'FontSize', 14); grid minor; grid on
                        set( ax(5), 'Tag', 'ax5');

                        % VarX vs Iteration

                    case 'iter'
                        % Simulated Image

                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        set( f, figProps{:});
    %                     subplot(232)
                        axes( findobj('Tag', 'ax2') );
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = findobj( get( gca,'Children'), 'Tag', 'imgSim');
                        set( img, 'CData', max( imageSim, [], 3) );
                        
                        % Residual Image
    %                     subplot(233)
                        axes( findobj('Tag', 'ax3') );
                        imgRes = max( abs( ImageOrg - imageSim), [], 3); 
                        img = findobj( get( gca,'Children'), 'Tag', 'imgRes');
                        set( img, 'CData', max( imgRes, [], 3) );

                        % Graphical Features
    %                     subplot(234)
                        axes( findobj('Tag', 'ax4') );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))

                        % Update residual plot
    %                     subplot(235)
                        axes( findobj('Tag', 'ax5') );
                        plotRes = findobj( get( gca,'Children'), 'Tag', 'plotRes');
                        X = [ get( plotRes, 'Xdata') optim.iteration]; Y = [ get( plotRes, 'Ydata') optim.resnorm ];
                        set( plotRes, 'Xdata', X, 'Ydata', Y);
                        set( get( gca, 'Title'), 'String', sprintf('Best Resnorm : %g',optim.resnorm) );

                        drawnow
                        pause(0.5)

                    case 'done'
                        % No clean up tasks required for this plot function.        
                        stop = true;
                end    

                % Make a folder and save these figures
                sName = [fitInfo.saveDirectory, filesep, sName];
                export_fig( sName, '-png', '-nocrop', '-a1') 
                %  }}}

            elseif graphicsVerbose == 0

                % do all this!
                % what should i save?
                % init:
                %   time
                %   channel
                %   fitScope
                %   featureNumber
                %   featureObj
                %
                % vector labels (1st iteration)
                % vector values (all iterations) (make an array)
                % 

            end


        end
        % }}}
        % save_midfit {{{
        function stop = save_midfit( vec, optim, state, fitInfo )
            % Plots during Lsqnonlin fitting
            % Inputs :  
            %   OPTIMVALUES: properties are funccount, iteration
            %   STATE: string representing current state of optimization engine. Values are init, iter, done 
            %   STOP: A boolean to stop the algorithm.
            %   All other inputs available to error function are also available. Currently: ImageOrg( original image), fitInfo and parameters

            stop = false;

            ImageOrg = fitInfo.image;
            nX = fitInfo.numVoxels(2); nY = fitInfo.numVoxels(1); nZ = fitInfo.numVoxels(3); 
            dim = length( size(ImageOrg) );
            image2D = max( ImageOrg, [], 3); intMin = min( ImageOrg(:) ); intMax = max( ImageOrg(:) );
            
            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

            time = fitInfo.time;
            channel = fitInfo.channel;
            fitScope = fitInfo.fitScope;
            if strcmp( fitScope, 'local')
                fName = sprintf('C%d_T%d_%s_feat%d_iter%d', channel, time, fitScope, fitInfo.featureIndex, optim.iteration );
                sName = sprintf('%s_feat%d_iter%d',  fitScope, fitInfo.featureIndex, optim.iteration );
            elseif strcmp( fitScope, 'global')
                fName = sprintf('C%d_T%d_%s_iter%d', channel, time, fitScope, optim.iteration );
                sName = sprintf('%s_iter%d', fitScope, optim.iteration );
            elseif strcmp( fitScope, 'globum_add')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            elseif strcmp( fitScope, 'globum_remove')
                fName = sprintf('C%d_T%d_globum%d-%d_iter%d', channel, time, fitInfo.Nold, fitInfo.Nnew, optim.iteration );
                sName = sprintf('globum%d-%d_iter%d', fitInfo.Nold, fitInfo.Nnew, optim.iteration );
            end
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            graphicsVerbose = 1;
            if graphicsVerbose == 1
                %  Interactive Plotting {{{
                switch state
                    case 'init'

                        % Make Figure
                        figTitle = [];
                        figName = [];
                        figPath = [];
                        f = figure( figProps{:} );
                        set( f, 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope) );
                        ax = tight_subplot(2,3, 0.05);
                        delete( ax(6) );
                        drawnow
                        pause(1)                        

                        % Original Image
                        axes( ax(1) ); %                     subplot(231)
                        img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');
                        set( ax(1), 'Tag', 'ax1');

                        % Simulated Image
    %                     subplot(232)
                        axes( ax(2) ); 
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = imagesc( max(imageSim, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated'); set( img, 'Tag', 'imgSim');
                        set( ax(2), 'Tag', 'ax2');

                        % Residual Image
    %                     subplot(233)
                        axes( ax(3) ); 
                        img = imagesc( max( abs( ImageOrg - imageSim), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Residual'); set( img, 'Tag', 'imgRes');
                        set( ax(3), 'Tag', 'ax3');

                        % Graphical Features
    %                     subplot(234)
                        axes( ax(4) ); 
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))
                        set( ax(4), 'Tag', 'ax4');

                        % Residual vs Iteration
    %                     subplot(235)
                        axes( ax(5) ); 
                        xtickformat( '%.2f'); ytickformat( '%d')
                        plotRes = plot( optim.iteration, optim.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set( plotRes,'Tag','plotRes');
                        xlabel('Iteration','interp','none'); ylabel('Resnorm','interp','none');
                        title( sprintf('Best Resnorm : %g', optim.resnorm ),'interp','none');
                        set( gca, 'FontSize', 14); grid minor; grid on
                        set( ax(5), 'Tag', 'ax5');

                        % VarX vs Iteration

                    case 'iter'
                        % Simulated Image

                        f = findobj( 'Tag', sprintf('fig %d_%d_%s', channel, time, fitScope));
                        set( f, figProps{:});
    %                     subplot(232)
                        axes( findobj('Tag', 'ax2') );
                        imageSim = Cell.simulateCellWithVec( vec, fitInfo);
                        img = findobj( get( gca,'Children'), 'Tag', 'imgSim');
                        set( img, 'CData', max( imageSim, [], 3) );
                        
                        % Residual Image
    %                     subplot(233)
                        axes( findobj('Tag', 'ax3') );
                        imgRes = max( abs( ImageOrg - imageSim), [], 3); 
                        img = findobj( get( gca,'Children'), 'Tag', 'imgRes');
                        set( img, 'CData', max( imgRes, [], 3) );

                        % Graphical Features
    %                     subplot(234)
                        axes( findobj('Tag', 'ax4') );
                        imagesc( image2D); eval(imageSets); hold on;
                        fitInfo.featureCurrent.displayFeature(gca);
                        set( get(gca, 'title'), 'String', sprintf('N = %d',fitInfo.featureMain.getSubFeatureNumber()))

                        % Update residual plot
    %                     subplot(235)
                        axes( findobj('Tag', 'ax5') );
                        plotRes = findobj( get( gca,'Children'), 'Tag', 'plotRes');
                        X = [ get( plotRes, 'Xdata') optim.iteration]; Y = [ get( plotRes, 'Ydata') optim.resnorm ];
                        set( plotRes, 'Xdata', X, 'Ydata', Y);
                        set( get( gca, 'Title'), 'String', sprintf('Best Resnorm : %g',optim.resnorm) );

                        drawnow
                        pause(0.5)

                    case 'done'
                        % No clean up tasks required for this plot function.        
                        stop = true;
                end    

                % Make a folder and save these figures
                sName = [fitInfo.saveDirectory, filesep, sName];
                export_fig( sName, '-png', '-nocrop', '-a1') 
                %  }}}

            elseif graphicsVerbose == 0

                % do all this!
                % what should i save?
                % init:
                %   time
                %   channel
                %   fitScope
                %   featureNumber
                %   featureObj
                %
                % vector labels (1st iteration)
                % vector values (all iterations) (make an array)
                % 

            end


        end
        % }}}
        % displayFinalFit {{{
        function h = displayFinalFit( Image2Fit, mainFeature, fitInfo)
            % Display Final Features and Fitting Results

            % Check if mainFeature has been dit. If not, then return out of method
            if isempty( mainFeature.featureList)
                fprintf('                       Skipping display for %s\n', mainFeature.type);  
                return
            end

            nX = size( Image2Fit, 2); nY = size( Image2Fit, 1); nZ = size( Image2Fit, 3); dim = length( size(Image2Fit) );
            image2D = max( Image2Fit, [], 3); intMin = min( Image2Fit(:) ); intMax = max( Image2Fit(:) );
            
            fName = sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope );
            figProps = {'NumberTitle', 'off', 'Name', fName, 'Position', [1 1 1280 720]};

            imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';
            h = figure( figProps{:} );
            ax = tight_subplot(2, 3, 0.05);
            drawnow; pause(3)

            % Original Image
%             subplot(231)
            axes( ax(1) );
            img = imagesc( image2D ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Original'); set( img, 'Tag', 'imgOrg');

            % Initial Simulated Image
%             subplot(232)
            axes( ax(2) );
            imageSimI = Cell.simulateCellWithVec( fitInfo.fitVecs.vec, fitInfo);
            img = imagesc( max(imageSimI, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated Init'); set( img, 'Tag', 'imgSimI');

            % Final Simulated Image
%             subplot(233)
            axes( ax(3) );
            imageSimF = Cell.simulateCellWithVec( fitInfo.fitResults.vfit, fitInfo);
            img = imagesc( max(imageSimF, [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', 'Image Simulated Final'); set( img, 'Tag', 'imgSimF');

            % Features
%             subplot(234) 
            axes( ax(4) );
            img = imagesc( image2D ); eval( imageSets); set( img, 'Tag', 'imgFeatures'); hold on;
            mainFeature.displayFeature(gca);
            set( get(gca, 'title'), 'String', sprintf('Features : N = %d',fitInfo.featureMain.getSubFeatureNumber()))

            % Initial Residual
%             subplot(235)
            axes( ax(5) ); 
            img = imagesc( max( abs( Image2Fit - imageSimI), [], 3) ); eval( imageSets); 
            set( get(gca, 'title'), 'String', 'Image Residual Init 100%'); set( img, 'Tag', 'imgResI');

            % Final Residual
%             subplot(236)
            axes( ax(6) );
            residFrac = 100*sum( abs( Image2Fit(:) - imageSimF(:) ) ) / sum( abs( Image2Fit(:) - imageSimI(:) ) );
            img = imagesc( max( abs( Image2Fit - imageSimF), [], 3) ); eval( imageSets); set( get(gca, 'title'), 'String', sprintf('Image Residual Final %.1f%%', residFrac) ); set( img, 'Tag', 'imgResF');

            sName = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ) ];
            export_fig( sName, '-png', '-nocrop', '-a1')

        end
        % }}}
        % dispImg {{{        
        function dispImg( varargin)
                
            if length(varargin) == 1 && ( isa( varargin{1}, 'double') || isa( varargin{1}, 'logical') || isa( varargin{1}, 'uint8') || isa( varargin{1}, 'single') )
                % its a single image
                img = varargin{1};
                figure; imagesc(img); axis equal; 
                if size(img,3) == 1
                    colormap gray
                end

                pos = get(gcf, 'position');
                set(gcf, 'pos', [-10000 pos(2:end)]);
                set(gcf, 'WindowState', 'maximized');
                set(gca, 'xlim', [1 size(img, 2)], 'ylim', [1 size(img, 1)], 'xtick',[], 'ytick',[] )
                
            elseif length( varargin) > 1
                % there are multiple images for comparison. The display order is
                % specified in the last argument [nrows ncols]
                
                if size( varargin{end}, 2) ~= 2
                    error('dispImg: the last argument must specify the order of plots [nRows nCols]')
                end
                
                nFig = length( varargin)-1;
                order = varargin{ end};
                nR = order(1); nC = order(2);
                img = varargin(1 : nFig);
                figure; 
                
                for jFig = 1 : nFig

                    subplot( nR, nC, jFig); imagesc( img{jFig} ); axis equal; 
                    if size( img{jFig}, 3) == 1
                        colormap gray
                    end
                    set(gca, 'xlim', [1 size(img{jFig}, 2)], 'ylim', [1 size(img{jFig}, 1)], 'xtick',[], 'ytick',[] )
                end
                pos = get(gcf, 'position');
                set(gcf, 'pos', [-10000 pos(2:end)]);
                set(gcf, 'WindowState', 'maximized');
                    
            else
                error('Oops, something went wrong!')
            end

        end
        % }}}

        % }}}
    
        % saveFinalFit {{{
        function h = saveFinalFit( Image2Fit, mainFeature, fitInfo)
            % Save data for this particular fit

            % Book-keeping
            channel = uint8(fitInfo.channel);
            time = uint16(fitInfo.time);
            fitScope = fitInfo.fitScope;
%             Image2Fit = uint16 (Image2Fit);

            % Images Simulated
            if ~isempty( mainFeature.featureList)
                imageSimI = uint16( Cell.simulateCellWithVec( fitInfo.fitVecs.vec, fitInfo) );
                imageSimF = uint16( Cell.simulateCellWithVec( fitInfo.fitResults.vfit, fitInfo) );
            else
                imageSimI = [];
                imageSimF = [];
            end

            % Main feature
            featureMainStruct = fitInfo.featureMain.saveAsStruct();

            f_data = [ fitInfo.saveDirectory, filesep,  sprintf('C%d_T%d_%s', fitInfo.channel, fitInfo.time, fitInfo.fitScope ), '.mat' ];
            save( f_data, 'channel', 'time', 'fitScope', 'imageSimI', 'imageSimF', 'featureMainStruct', '-v7');

        end
        % }}}
        
    end

    methods ( Static = true , Access = protected )

        % parseSettings {{{
        function settings = parseSettings( settings)
            
            flags.doFit = settings.CFGinfo.runFit;
            flags.graphRealTime = settings.CFGinfo.runRealTimeGraphics;
            flags.graphPostFit = settings.CFGinfo.runPostFitGraphics;
            flags.track = settings.CFGinfo.runTracking;
            flags.analyze = settings.CFGinfo.runFit;
            flags.debug = settings.CFGinfo.runDebug; % enables all graphics
            flags.doFitLocal = settings.CFGinfo.fitLocal;
            flags.doFitFeatureNumber = settings.CFGinfo.fitFeatureNumber;
            flags.fitExploreSpeed = settings.CFGinfo.fitExploreSpeed;
            if isfield( settings, 'fitSpindleOnly')
                flags.fitSpindleOnly = 1;
            else
                flags.fitSpindleOnly = 0;
            end
            if isfield( settings, 'skipOptimizeNumber')
                flags.skipOptimizeNumber = settings.skipOptimizeNumber;
            else
                flags.skipOptimizeNumber = 0;
            end
            if isfield( settings, 'removeSpindleSPB')
                flags.removeSpindleSPB = settings.removeSpindleSPB;
            else
                flags.removeSpindleSPB = 0;
            end

            settings.paths.saveParent = settings.CFGinfo.saveParent;
            settings.paths.run = settings.CFGinfo.runPath;
            settings.paths.movie = settings.cellinfo.moviePath;
            settings.flags = flags;

        end
        % }}}

    end

    methods ( Abstract = true )
    
        obj = findFeatures( obj, image, time, channelTrue, channelIdx)

    end

end
