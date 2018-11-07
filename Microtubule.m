classdef Microtubule
    % This defines a Microtubule    
    properties
        % properties {{{
        mtoc % nucleation point
        orientationInit % approx orientation (in radians)
        source2D % 2D source image
        source % source image containing the microtubule
        mask % mask for microtubule
        mask2D % mask in 2D
        helperImage % helper image containing just the microtubule intensity
        display % display features specific to this microtubule (color, linewidth, markersize, markerstyle, etc...)
        estimatedCoords % the estimated points through which the microtubule passes
        polyOrder % polynomial order( 2=quadratic, 3=cubic)
        estimatedCoef % estimated polynomial coefficients
        amplitude % gaussian amplitude
        std % gaussian standard deviation
        background % background intensity
        fitCoef % fitted polynomial coefficients
        id % identification
        fitProps % properties related to fitting
        dead % is the microtubule alive or dead. dead microtubules are discarded and data is lost.
        dim % dimensionality of the microtubule
        % }}}
    end
    
    % Main Methods
    methods
        
        % Initialization
        % Microtubule {{{
        function obj = Microtubule( sourceImage, mask, displayImage, coords, orient)
            % To create a microtubule, all you need is a source image and a
            % nucleation point (x and y)
            obj.dim = size( coords, 1);
            obj.mtoc = coords(:,1);
            obj.estimatedCoords.x = coords(1,:);
            obj.estimatedCoords.y = coords(2,:);
            if obj.dim==3; obj.estimatedCoords.z = coords(3,:); end
            obj.orientationInit = orient;
            obj.source = sourceImage;
            obj.source2D = max( sourceImage, [], 3);
            obj.mask = mask;
            obj.mask2D = max( mask, [], 3);
            obj.display.MarkerSize = 12;
            obj.display.LineWidth = 4;
            obj.dead = 0;
            obj.display.image = displayImage;
        end
        % }}} 
        
        % Estimation
        % measureIntensity {{{
        function intensities = measureIntensity( obj, type)
            % used to measure intensity while using estimated coefficients
            t = linspace(0,1);
            
            if strcmp(type, 'estimate')
                px = obj.estimatedCoef(1,:);
                py = obj.estimatedCoef(2,:);
                
            elseif strcmp(type, 'fit')
                px = obj.fitCoef(1,:);
                py = obj.fitCoef(2,:);
                
            end
            x = polyval(px, t);
            y = polyval(py, t);
            
            coord = round( [x; y]);
            % keep nique coordinates
            idxRm = [];
            for jCrd = 1 : size(coord, 2)-1
                if ( coord(1, jCrd) == coord(1, jCrd+1) ) && ( coord(2, jCrd) == coord(2, jCrd+1) )
                    idxRm = [idxRm, jCrd];
                end
            end

            coord(:, idxRm) = [];
            idx = sub2ind( size(obj.source2), coord(2, :), coord(1, :) );
            intensities = obj.source2( idx);
            
        end
        % }}}

        % EstimateMicrotubuleCurve {{{
        function obj = EstimateMicrotubuleCurve(obj, polyOrder)
            % fits a polynomial of order 'polyOrder' to obj.estimatedPoints
            
            obj.polyOrder = polyOrder;
            
            % Get data to fit
            x = obj.estimatedCoords.x;
            y = obj.estimatedCoords.y;
            
            t = linspace(0,1);
            % interpolate for better accuracy
            xi = interp1( linspace(0,1,length(x) ), x, t);
            yi = interp1( linspace(0,1,length(y) ), y, t);
            
            % fit polynomial of order polyOrder
            if obj.dim == 2
                obj.estimatedCoef = [ polyfit( t, xi, polyOrder); polyfit( t, yi, polyOrder) ];
            elseif obj.dim == 3
                z = obj.estimatedCoords.z;
                zi = interp1( linspace(0,1, length(z) ), z, t);
                obj.estimatedCoef = [ polyfit( t, xi, polyOrder); polyfit( t, yi, polyOrder) ; plotfit( t, zi, polyOrder) ];
            end

        end
        % }}}

        % EstimateGaussianParameters {{{        
        function obj = EstimateGaussianParameters( obj, params)

            obj.amplitude = mean( measureIntensity( obj, 'estimate') );
            obj.std = params.structInit.std;;
            obj.background = median( obj.source2( obj.mask2D) );

        end
        % }}}

        % Fitting
        % fitInitialization {{{
        function obj = fitInitialization(obj, params) 

            params.initialization.initVals(1:4) = [ obj.background, obj.amplitude, obj.stdev ];
            params.initialization.ub = params.initialization.initVals(1:4) + params.initialization.vecRange;
            params.initialization.lb = params.initialization.initVals(1:4) - params.initialization.vecRange;
            
            params.initialization = initializeParameterStructures( obj, params.initialization);
            [vec, vecUb, vecLb] = convertStruct2Vec_MT( params.initialization, params.fixMTOC);

            % create optimization options and the problem for lsqnonlin
            opts = optimoptions( @lsqnonlin, 'MaxFunEvals', params.maxFunEvals, 'OptimalityTolerance', params.optTol, 'MaxIter', params.maxIter, ...
            'TolFun', params.tolFun, 'FiniteDifferenceStepSize', params.finDifStepSize, 'FiniteDifferenceType', 'central', 'StepTolerance', params.stepTol, ...
            'display', params.display, 'OutputFcn', @plotOptimStatus);
    
            % parameters for fitting function
            FixedParams.PolyOrder = params.polyOrder;
            FixedParams.StartingPointFixed = params.fixMTOC;
            FixedParams.XCoefEnd = obj.estimatedCoef(1,end);
            FixedParams.YCoefEnd = obj.estimatedCoef(2,end);
            FixedParams.Mask2D = obj.mask2D;
            FixedParams.Mask = obj.mask;
            FixedParams.dim = params.fitDim;
            FixedParams.numberOfMicrotubules = 1; % local fitting

            if params.fitDim == 2
                FixedParams.Mask = 1.0 * logical(obj.source2);
                ImageFit = obj.source2.*obj.mask2D;
            elseif params.fitDim == 3
                FixedParams.Mask = 1.0 * logical( obj.source);
                ImageFit = obj.source.*obj.mask;
            end

            % make the error function for lsqnonlin
            f = makeErrorFcn( ImageFit, FixedParams, params);

            % Crate Optimization Problem
            prob = createOptimProblem( 'lsqnonlin', 'objective', f, 'x0', vec, 'ub', vecUb, 'lb', vecLb, 'options', opts);
            
            % Save problem to object for easy access
            obj.fitProps.problem = prob;

            % initializeParameterStructures {{{
            function params = initializeParameterStructures( obj, params) 
                % we will allow the curve length to be altered between 0.5 and 2 (and use the coef rang eas upper and lower bounds (which will give it some freedom)

                % we can do this easily by parameterizing our curve values to be representative of values from t=0 to t = t_max (where t_max is 1 by default). Then we just allow t_max to lie between 0.5 and 2.5 and we fit curves again through it using the default t-range and we end up with a range of coeficients.
                fitdim = params.structInit.dim;
                polyOrder = params.structInit.polyOrder;
                lenMin = params.lengthMin;
                lenMax = params.lengthMax;
                sigKeep = 0.5;
                tRange = linspace( lenMin, lenMax, 100);
             
                coefs = obj.estimatedCoef;
                tDef = linspace( 0, 1);
                xcfs = []; ycfs = []; zcfs = [];
                for jj = 1 : length(tRange)
                    tNew = linspace( 0, tRange(jj) );
                    xNew = polyval( coefs( 1, :), tNew);
                    yNew = polyval( coefs( 2, :), tNew);
                    if fitdim==3, zNew = polyval( coefs( 3, :), tNew); end

                    xcfs = [xcfs ; polyfit( tDef, xNew, polyOrder) ];
                    ycfs = [ycfs ; polyfit( tDef, yNew, polyOrder) ];
                    if fitdim==3, zcfs = [zcfs ; polyfit( tDef, zNew, polyOrder) ]; end
                end

                % find the range of these coefs ( model with a gaussian and keep up to a certain number of standard deviations)
                sigX = std( xcfs, 0, 1);
                sigY = std( ycfs, 0, 1);
                if fitdim==3, sigZ = std( zcfs, 0, 1); end

                % We will keep half a standard deviation above the max coeff val and half a std below the min coeff value
                ubX = max( xcfs, [], 1) + sigKeep * sigX;
                ubY = max( ycfs, [], 1) + sigKeep * sigY;
                lbX = min( xcfs, [], 1) - sigKeep * sigX;
                lbY = min( ycfs, [], 1) - sigKeep * sigY;
                if fitdim==3, 
                    ubZ = max( zcfs, [], 1) + sigKeep * sigZ;
                    lbZ = min( zcfs, [], 1) - sigKeep * sigZ;
                end 

                params.structInit.coefX{1} = coefs(1,:);
                params.structInit.coefY{1} = coefs(2,;);
                if fitdim==3, bounds.structInit.coefZ{1} = coefs(3,:); end
                
                params.structUb.coefX{1} = ubX;
                params.structUb.coefY{1} = ubY;
                params.structLb.coefX{1} = lbX;
                params.structLb.coefY{1} = lbY;
                if fitdim==3, params.structUb.coefZ{1} = ubZ; params.structLb.coefZ{1} = lbZ; end
                
                params.structInit.background = obj.background;
                params.structUb.background = obj.background + params.structUb.background;
                params.structLb.background = obj.background + params.structLb.background;

                params.structInit.amplitude= obj.amplitude;
                params.structUb.amplitude = obj.amplitude + params.structUb.amplitude;
                params.structLb.amplitude= obj.amplitude + params.structLb.amplitude;
                
            end
            % }}}

            % makeErrorFcn {{{    
            function errVal = makeErrorFcn( ImageFit, FixedParams, params)

                errVal = @ErrorFcn;

                function err = ErrorFcn( p)
                        
                    err = feval( params.fitfunc, p, ImageFit, FixedParams) - ImageFit;
                    % err = GaussianCurveMaker2D( p, Image2D, FixedParams) - Image2D;
                    err = err(:);

                end

            end

            % }}}

                % plotOptimStatus{{{
                function stop = plotOptimStatus( x, optimValues, state)
                % plotOptimStatus: updates a plot at every optim iteration showing optimization in real-time
                %   OPTIMVALUES: Information after the current local solver call.
                %          funccount: number of function evaluations
                %          iteration: iteration number
                %   STATE: Current state in which plot function is called. 
                %          Possible values are:
                %             init: initialization state 
                %             iter: iteration state 
                %             done: final state
                %   STOP: A boolean to stop the algorithm.
                %   Copyright 2014 The MathWorks, Inc.

                % Initialize stop boolean to false.
                stop = false;
                numPixX = size(FixedParams.Mask,1);
                numPixY = size(FixedParams.Mask,2);

                if length( size(ImageFit) ) == 3
                    Image2D = max( ImageFit, [], 3);
                elseif length( size(ImageFit) ) == 2
                    Image2D = ImageFit;
                end

                switch state
                    case 'init'
                        % init {{{
                        t = linspace(0,1); % paramateric variable
                        startIdx = 5;
                        gap = FixedParams.PolyOrder;
                        if FixedParams.StartingPointFixed
                            xcurr = polyval( [ x( startIdx : startIdx+ gap-1), FixedParams.XCoefEnd], t); % current curve coordinates
                            ycurr = polyval( [ x( startIdx+ gap : end), FixedParams.YCoefEnd], t); % current curve coordinates
                        else
                            xcurr = polyval( x( startIdx : startIdx+ gap-1), t); % current curve coordinates
                            ycurr = polyval( x( startIdx+ gap : end), t); % current curve coordinates
                        end

                        figure('Name', ['lsqnonlin_', num2str(obj.id)], 'NumberTitle', 'off');
                        posOld = get(gcf, 'Position');
                        set(gcf, 'position', [-10000 posOld(2:end)], 'WindowState', 'maximized');
                        subplot(121)
                        plotBest = plot(optimValues.iteration,optimValues.resnorm, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set(plotBest,'Tag','psoplotbestf');
                        xlabel('Iteration','interp','none');
                        ylabel('Function value','interp','none')
                        title(sprintf('Best Function Value: %g',optimValues.resnorm ),'interp','none');
                        set(gca, 'FontSize', 14)
                        grid minor; grid on

                        subplot(122)
                        imagesc([Image2D]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.5], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow
                        % }}}
                    case 'iter'
                        %  iter {{{
                        subplot(121)
                        plotBest = findobj(get(gca,'Children'),'Tag','psoplotbestf');
                        newX = [get(plotBest,'Xdata') optimValues.iteration];
                        newY = [get(plotBest,'Ydata') optimValues.resnorm ];
                        set(plotBest,'Xdata',newX, 'Ydata',newY);
                        set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimValues.resnorm) );
                        grid minor; grid on

                        t = linspace(0,1); % paramateric variable
                        startIdx = 5;
                        gap = FixedParams.PolyOrder;
                        if FixedParams.StartingPointFixed
                            xcurr = polyval( [ x( startIdx : startIdx+ gap-1), FixedParams.XCoefEnd], t); % current curve coordinates
                            ycurr = polyval( [ x( startIdx+ gap : end), FixedParams.YCoefEnd], t); % current curve coordinates
                        else
                            xcurr = polyval( x( startIdx : startIdx+ gap-1), t); % current curve coordinates
                            ycurr = polyval( x( startIdx+ gap : end), t); % current curve coordinates
                        end

                        subplot(122)
                        imagesc([Image2D]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow
                        % }}}
                    case 'done'
                        % No clean up tasks required for this plot function.        
                end    

                end
                % }}}
           
            % outputFcn {{{
            function stop = outputFcn( x, optimValues, state)
            % outputFcn : code that is run after every successful iteration of the optimzation solver
            % x : current state vector
            % 
            stop = false;

            switch state
                case 'init'
                    labels = cell( length( x) );
                    labels{1} = 'bkg';
                    labels{2} = 'amp';
                    labels{3} = 'sigX';
                    labels{4} = 'sigY';
                    labels{5} = 'polyCoef';
                    % Create labels for the parameter vectors and save the values in x to a csv file
                    print('output function')
                case 'iter'
                    print('output function') 
                case 'done'
                % Cleanup of plots, guis, or final plot
            end
            % }}}
            
            end
        
        end
        % }}}
        
        % fitMicrotubule{{{
        function obj = fitMicrotubule( obj, params)
            
            
            [vfit,resnorm,residual,exitflag,~,~,~] = lsqnonlin( obj.fitProps.problem); 

            pStructFin.bkg = vfit(1);
            pStructFin.amp = vfit(2);
            pStructFin.sigma = vfit(3);

            startIdx = 5;
            if params.fixMTOC
                pStructFin.XCoef = [ vfit( startIdx: startIdx+obj.polyOrder-1), obj.estimatedCoef(1,end)];
                pStructFin.YCoef = [ vfit( startIdx+obj.polyOrder :end), obj.estimatedCoef(2,end)];
            else
                pStructFin.XCoef = [ vfit( startIdx: startIdx+obj.polyOrder-1) ]; 
                pStructFin.YCoef = [ vfit( startIdx+obj.polyOrder :end) ];
            end
            
            obj.fitProps.vfit = vfit;
            obj.fitProps.structFit = pStructFin
            obj.fitCoef = [pStructFin.XCoef ; pStructFin.yCoef];
            
        end
        % }}}

        % Plotting
        % PlotMicrotubuleCurve {{{
        function PlotMicrotubuleCurve(obj, type)
            % plots the polynomial coming from 'type' = 'estimate' or 'fit'
            
            t = linspace(0,1);
            
            if strcmp(type, 'estimate')
                px = obj.estimatedCoef(1,:);
                py = obj.estimatedCoef(2,:);
                
                figName = sprintf( 'mt_%d_poly_%d_estimate', obj.id, obj.polyOrder);
                figTitle = sprintf( 'Microtubule %d poly-%d estimate,', obj.id, obj.polyOrder);
            elseif strcmp(type, 'fit')
                px = obj.fitCoef(1,:);
                py = obj.fitCoef(2,:);
                
                figName = sprintf( 'mt_%d_poly_%d_fit', obj.id, obj.polyOrder);
                figTitle = sprintf( 'Microtubule %d poly-%d fit,', obj.id, obj.polyOrder);
            end
            
            figure('Name', figName, 'NumberTitle', 'off')
            pos = get(gcf, 'position');
            set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
            imagesc( obj.display.image); axis equal; colormap gray;
            set(gca, 'xlim', [1 size( obj.source,1)], 'ylim', [1 size( obj.source,2)]); hold on;

            polyOrder = obj.polyOrder;
            x = polyval(px, t);
            y = polyval(py, t);

            plot( x, y, 'r-', 'LineWidth', obj.display.LineWidth ); 
            plot( x(1), y(1), 'r*', 'MarkerSize', obj.display.MarkerSize+3, 'LineWidth', obj.display.LineWidth/2);
            LH(1) = plot(nan, nan, 'r-', 'LineWidth', obj.display.LineWidth); L{1} = ['MT ', num2str(obj.id)];
            LH(2) = plot(nan, nan, 'r*', 'MarkerSize', obj.display.MarkerSize); L{2} = 'MTOC';
            legend(LH, L);
            set(gca, 'FontSize', 15);
            title( figTitle);
            hold off; 
            set(gcf, 'WindowState', 'maximized');
            
        end
        % }}}

        % plotEstimatedPointNext {{{
        function plotEstimatedPointNext( obj, orientation, stepSize, visibility, fieldOfVision, iter)
            % Plot next estimated point stored in obj.estimatePoints and
            % overlays it on obj.source. Also plots visibility circle and
            % field of vision arc
            if nargin == 6
                figName = sprintf( 'mt_%d_estimate_iter_%d', obj.id, iter);
                figTitle = sprintf( 'Microtubule %d estimate : iteration %d, param=(%.1f,%.1f),', obj.id, iter, stepSize, visibility);
            else
                figName = sprintf( 'mt_%d_estimate', obj.id);
                figTitle = sprintf( 'Microtubule %d estimate, param=(%.1f,%.1f),', obj.id, stepSize, visibility);
            end
            
            figure('Name', figName, 'NumberTitle', 'off')
            pos = get(gcf, 'position');
            set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
            imagesc( obj.display.image); axis equal; colormap gray;
            set(gca, 'xlim', [1 size( obj.source,1)], 'ylim', [1 size( obj.source,2)]); hold on;
            
            % create visibility circle and plot field of Vision in a
            % different color
            x0 = obj.estimatedPoints(1,end-1);
            y0 = obj.estimatedPoints(2,end-1);
            
            x1 = obj.estimatedPoints(1,end);
            y1 = obj.estimatedPoints(2,end);
            
            % create boundary for visibility circle
            th = 0:pi/50:2*pi;
            xCir = visibility * cos(th) + x0;
            yCir = visibility * sin(th) + y0;
            
            % create boundary for field of vision arc
            angLimits = orientation + deg2rad( [1,-1]*fieldOfVision );
            th2 = angLimits(2) : pi/50: angLimits(1);
            xCir2 = visibility * cos(th2) + x0;
            yCir2 = visibility * -sin(th2) + y0;
            
            % create fill region for field of vision
            xx = [linspace(x0, xCir2(1),30), xCir2, linspace(xCir2(end), x0,30)]';
            yy = [linspace(y0, yCir2(1),30), yCir2, linspace(yCir2(end), y0,30)]';
            
            % plot vision torchlight
            fill( xx, yy, 'y');
            
            % plot all estimated points with links in between them
            plot( obj.estimatedPoints(1,1:end-1), obj.estimatedPoints(2,1:end-1), 'r-', 'LineWidth', obj.display.LineWidth, ...
                'Marker', '.', 'MarkerSize', obj.display.MarkerSize)
            
            % plot current search point
            plot( obj.estimatedPoints(1,end-1), obj.estimatedPoints(2,end-1), 'rx', 'LineWidth', obj.display.LineWidth, ...
                'MarkerSize', obj.display.MarkerSize)
            
            % plot mtoc location
            plot( obj.estimatedPoints(1,1), obj.estimatedPoints(2,1), 'r*', 'LineWidth', obj.display.LineWidth-1, ...
                'MarkerSize', obj.display.MarkerSize+4)
            
            % plot visibility circle, field of vision arc
%             plot(xCir, yCir, 'b','LineWidth', obj.display.LineWidth/2); % plot visib circle
%             plot(xCir2, yCir2, 'g','LineWidth', obj.display.LineWidth/2); % plot field of vision arc
            
            % plot next possible link
            plot( obj.estimatedPoints(1,end), obj.estimatedPoints(2,end), 'ro', 'LineWidth', obj.display.LineWidth, ...
                'MarkerSize', obj.display.MarkerSize)
            
            % create legend entry
            LH(1) = plot(nan, nan, 'rx', 'MarkerSize', obj.display.MarkerSize); L{1} = 'Search Point';
            LH(2) = plot(nan, nan, 'ro', 'MarkerSize', obj.display.MarkerSize); L{2} = 'Next Link';
%             LH(3) = plot(nan, nan, 'b-', 'MarkerSize', obj.display.MarkerSize); L{3} = 'Forbidden Angles';
%             LH(4) = plot(nan, nan, 'g-', 'MarkerSize', obj.display.MarkerSize); L{4} = 'Allowed Angles';
            LH(3) = plot(nan, nan, 'r-', 'MarkerSize', obj.display.MarkerSize); L{3} = 'Estimated MT';
            LH(4) = plot(nan, nan, 'r*', 'MarkerSize', obj.display.MarkerSize); L{4} = 'MTOC';
            legend(LH, L);
            set(gca, 'FontSize', 15);
            title( figTitle);
            hold off
            
        end
        % }}}

        % plotEstimatedPointsAll {{{
        function plotEstimatedPointsAll( obj)
            % Plots all estimated points stored in obj.estimatePoints and
            % overlays it on obj.source

            figName = sprintf( 'mt_%d_estimate_points', obj.id);
            figTitle = sprintf( 'Microtubule %d estimate', obj.id);
            
            figure('Name', figName, 'NumberTitle', 'off')
            pos = get(gcf, 'position');
            set(gcf, 'pos', [-10000 pos(2:end)], 'WindowState', 'maximized');
            img = imagesc( obj.display.image); axis equal; colormap gray;
            set(gca, 'xlim', [1 size( obj.source,1)], 'ylim', [1 size( obj.source,2)]); hold on;
            
            x0 = obj.estimatedPoints(1,end);
            y0 = obj.estimatedPoints(2,end);
            
            % plot all estimated points with links in between them
            plot( obj.estimatedPoints(1,:), obj.estimatedPoints(2,:), 'r-', 'LineWidth', obj.display.LineWidth, ...
                'Marker', '.', 'MarkerSize', obj.display.MarkerSize)
            
            % plot mtoc location
            plot( obj.estimatedPoints(1,1), obj.estimatedPoints(2,1), 'r*', 'LineWidth', obj.display.LineWidth, ...
                'MarkerSize', obj.display.MarkerSize)
            
            % create legend entry
            LH(1) = plot(nan, nan, 'r-', 'MarkerSize', obj.display.MarkerSize); L{1} = 'Estimated MT';
            LH(2) = plot(nan, nan, 'r*', 'MarkerSize', obj.display.MarkerSize); L{2} = 'MTOC';
            legend(LH, L);
            set(gca, 'FontSize', 15);
            title( figTitle);
            hold off
            
        end
        % }}} 
        
        function outputArg = trackVariableChange(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

    end
end

