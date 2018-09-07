classdef Microtubule
    % This defines a Microtubule    
    properties
        mtoc % nucleation point
        orientationInit % approx orientation (in radians)
        source % source image containing the microtubule
        helperImage % helper image containing just the microtubule intensity
        display % display features specific to this microtubule (color, linewidth, markersize, markerstyle, etc...)
        estimatedPoints % the estimated points through which the microtubule passes
        polyOrder % polynomial order( 2=quadratic, 3=cubic)
        estimatedCoef % estimated polynomial coefficients
        amplitude % gaussian amplitude
        stdev % gaussian standard deviation
        background % background intensity
        fitCoef % fitted polynomial coefficients
        id % identification
        fitProps % properties related to fitting
        dead % is the microtubule alive or dead. dead microtubules are discarded and data is lost.
    end
    
    % Main Methods
    methods
        
        % Initialization
        % Microtubule {{{
        function obj = Microtubule( sourceImage, nucleateX, nucleateY, orient)
            % To create a microtubule, all you need is a source image and a
            % nucleation point (x and y)
            obj.mtoc = [ nucleateX, nucleateY];
            obj.orientationInit = orient;
            obj.source = sourceImage;
            obj.display.MarkerSize = 12;
            obj.display.LineWidth = 4;
            obj.dead = 0;
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
            idx = sub2ind( size(obj.source), coord(2, :), coord(1, :) );
            intensities = obj.source( idx);
            
        end
        % }}}

        % EstimateMicrotubuleCurve {{{
        function obj = EstimateMicrotubuleCurve(obj, polyOrder)
            % fits a polynomial of order 'polyOrder' to obj.estimatedPoints
            
            obj.polyOrder = polyOrder;
            
            % Get data to fit
            x = obj.estimatedPoints(1,:);
            y = obj.estimatedPoints(2,:);
            t = linspace(0,1);
            
            % interpolate for better accuracy
            xi = interp1( linspace(0,1,length(x) ), x, t);
            yi = interp1( linspace(0,1,length(y) ), y, t);
            
            % fit polynomial of order polyOrder
            obj.estimatedCoef = [ polyfit( t, xi, polyOrder); polyfit( t, yi, polyOrder) ];
            
        end
        % }}}

        % EstimateMicrotubulePoints {{{
        function obj = EstimateMicrotubulePoints(obj, stepSize, visibility, fieldOfVision)
            % EstimateMicrotubulePoints: estimates points along
            % microtubule.
            % We take the following approach:
            % 1) Nucleate at the MTOC. We find the angle that gives the maximum in a
            %    radially integrated signal, up to some max radius( which we will call
            %    the visibility)
            % 3) We'll propagate from the MTOC to a new point along the optimum angle,
            %    which a defined distance away from the MTOC (this we will call the
            %    step size). Throughout this process, we might look for optimum angles
            %    within a certain range related to the orientation of the Tube (this we
            %    will refer to as the field of vision).
            % 4) We will continue propagating until we can't find an optimum angle
            %    anymore (all our surroudings appear the same. This will be based on
            %    some optimality condition. At this point we will change our step size
            %    to half its original value and see if we can propogate a smaller
            %    distance forward.
            % 5) This will conclude our initial estimation phase.

            
            if isempty( obj.helperImage)
                helperImg = obj.source;
            else
                helperImg = obj.helperImage;
            end
            
            % clear out any previously stored estimated points
            obj.estimatedPoints = [];
            
            success = 1; iter = 1;
            max_mt_length = 80; % max MT length approx
            max_iter = round(max_mt_length/stepSize);
            orientationOld = obj.orientationInit;
            
            % iteratively propagate along microtubule
            while success && iter< max_iter
                
                [obj, success, orientationNext] = EstimateNextPoint( obj, orientationOld, stepSize, visibility, fieldOfVision, helperImg);
                if success
%                     plotEstimatedPointNext( obj, orientationNext, stepSize, visibility, fieldOfVision, iter)
%                     drawnow; pause(0.2)
                    orientationOld = orientationNext;

                end
                
                iter = iter+1;
                
            end
            
            % try continuing with half the stepSize
            success = 1; count=iter;
            
            while success && iter < count+5
                
                [obj, success, orientationNext] = EstimateNextPoint( obj, orientationOld, stepSize/2, visibility/2, fieldOfVision, helperImg);
            
                if success
%                     plotEstimatedPointNext( obj, orientationNext, stepSize, visibility, fieldOfVision, iter)
%                     drawnow; pause(0.2)
                    orientationOld = orientationNext;

                end
                iter = iter+1;
                
            end
            
            % if too short a guess, discard it
            if size(obj.estimatedPoints, 2) < 3
                obj.dead = 1;
            end
        end
        
        function [obj, success, phiFinal] = EstimateNextPoint( obj, orientation, stepSize, visibility, fieldOfVision, helperImage)
            % Given a starting point and some other parameters, this
            % function will find the next point that is connected by high
            % intensity.

            if ~isempty( obj.estimatedPoints)
                xPrev = obj.estimatedPoints(1,:);
                yPrev = obj.estimatedPoints(2,:);
            else
                xPrev = obj.mtoc(1);
                yPrev = obj.mtoc(2);
                obj.estimatedPoints(1,:) = xPrev;
                obj.estimatedPoints(2,:) = yPrev;
            end
            
            % This is the start point for propagation
            x0 = xPrev(end);
            y0 = yPrev(end);
            
            % We will create a finer image to allow us to use smaller step
            % sizes.
            [x, y] = meshgrid( 1 : size( helperImage, 1) );
            xVecNew = 1 : 0.25: size( helperImage, 1);
            [xi, yi] = meshgrid( xVecNew );
            imSub = interp2(x,y,helperImage,xi,yi,'linear');
            imSub2 = interp2(x,y,obj.display.image,xi,yi,'linear');

            clear x y xi yi
            
            angRange = deg2rad( fieldOfVision/2 );

            % Define the step sizes to use for the angular sweep in phi. We
            % start at 90 degrees to the proposed orientation of the tube
            PhiStep = deg2rad ( 2.5 ); % 2.5 degrees
            PhiVec = orientation - angRange : PhiStep : orientation + angRange;

            % Pre-allocate array for storing intensity data with varying phi.
            IntPhi = zeros( 1, length( PhiVec) );

            rp = visibility; % used for radial integration (allows for better angular determination)
            rmin = 1;
%             figure;
            for jPhi = 1 : length( PhiVec)
                % Extract the relevant angles for drawing a non-zero angle
                phi1 = PhiVec( jPhi) - PhiStep;
                phi2 = PhiVec( jPhi) + PhiStep;

                % Find the coordinates of the extreme points on the arc length
                % of this pizza slice.
                X1 = x0 + ( visibility * cos( [ phi1, phi2]) );
                Y1 = y0 + ( visibility * sin( [ phi1, phi2]) ); % -ve because matlab orientation is reversed( origin is at top left)

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
%                 imshowpair(imSub, imMask); pause(0.25)
                
                % Sum up values and store them
                IntPhi( jPhi) = sum( imMasked(:) ) / sum(imMask(:));
                
            end
            % We smooth the angular intensity to enable easy peak finding
            imSmooth = imgaussfilt( IntPhi, 3);
            
            if min(imSmooth) == max(imSmooth)
                success = 0; phiFinal = NaN;
                return
            end
            
            imSmooth = ( max(IntPhi)-min(IntPhi) ) * mat2gray( imSmooth ) + min( IntPhi);
            
            % Now we need to come up with a threshold for defining what is
            % and what isn't a peak. We'll do local intensity search for
            % this
            imMask = 0*helperImage;
            for jY = 1 : size(helperImage, 1)
                for jX = 1 : size(helperImage, 2)
                    if ( norm( [x0, y0] - [jX, jY] ) < rp ) && ( norm( [x0, y0] - [jX, jY] ) > rmin )
                        imMask( jY, jX) = 1;
                    end
                end
            end
            imMask = imMask .* helperImage;
            imPosVals = imMask( imMask > 0);
            T = multithresh( imPosVals, 1 );
            imSmooth( imSmooth < median(imPosVals) ) = median(imPosVals);
            if T < median(imPosVals); T = median(imPosVals); end
            
            
            if min(imSmooth) == max(imSmooth)
                success = 0; phiFinal = NaN;
                return
            end
            
            % Find Peaks Properties
            halfwidth = asin( 1.5/ stepSize);
            if ~isreal(halfwidth); halfwidth=asin(1); end;
            minProm = 0.01;
            minHeight = mean( [T, median(imPosVals)]);
            props = {'SortStr', 'descend', 'MinPeakProminence', minProm, 'MinPeakWidth', halfwidth, 'MinPeakHeight', minHeight};

            % Find the angle corresponding to the peaks in this region
            try;[pkInt, phi0Loc, ~, ~] = findpeaks( imSmooth, props{:} );catch
                stopH = 1;
            end
            phi0 = PhiVec( phi0Loc);
            
            % pick the brightest peak.
            [~, maxIdx] = max(pkInt);
            phiFinal = phi0( maxIdx);
            
            if ~ isempty(phiFinal)
                phiFinal = phiFinal';
                
                x1 = x0 + ( stepSize * cos( phiFinal ) );
                y1 = y0 + ( stepSize * sin( phiFinal ) );

                obj.estimatedPoints = [ obj.estimatedPoints(1,:), x1 ; obj.estimatedPoints(2,:), y1];
                success = 1;

            else
                success = 0;
            end

        end
        % }}}

        % EstimateGaussianParameters {{{        
        function obj = EstimateGaussianParameters( obj, params)

            obj.amplitude = measureIntensity( obj, 'estimate')
            if nargin < 2
                obj.stdev = params.sigma;
            else
                obj.stdev = [1.2, 1.2, 1.0];
            end
            obj.background = median( obj.source( obj.source > 0) );


        end
        % }}}

        % Fitting
        % prepareForFit {{{
        function obj = prepareForFit(obj, fitMethod, fixStartPoint)
            % assigns all required values to the microtubule and prepares
            % the upper, lower bound, options etc for the fitting method.
            % It takes in method = 'lsqnonlin' or 'particleswarm'
            
            % first, we need to measure the amplitude, background, and the
            % gaussian standard deviation
            obj.fitProps = [];
            
            intensities = measureIntensity( obj, 'estimate');
            try
                obj.amplitude = mean( intensities(end-4:end) );
            catch
                obj.amplitude = 0.25;
            end
            
            obj.background = median( obj.source( obj.source > 0) );
            obj.stdev = 1.2;
            
            if obj.amplitude <= obj.background
                warning( 'amplitude is less than background value for microtubule %d', obj.id)
            end
            
            % Define upper and lower limits on parameters
            ampRange = 0.02;
            stdevRange = 0.01;
            bkgrange = 0.01;
            
            AmpUb = obj.amplitude + ampRange; AmpLb = obj.amplitude - ampRange;
            SigUb = obj.stdev + stdevRange; SigLb = obj.stdev - stdevRange; 
            BkgUb = obj.background + bkgrange; BkgLb = obj.background - bkgrange;
            
            % ensure amplitude is always above background
            if AmpLb <= BkgUb
                warning( 'amplitude LB is less than background UB value for microtubule %d', obj.id)
            end
            
            if fixStartPoint == 1
                xP = obj.estimatedCoef(1, 1:end-1);
                yP = obj.estimatedCoef(2, 1:end-1);
            else
                xP = obj.estimatedCoef(1, :);
                yP = obj.estimatedCoef(2, :);
            end
            obj.fitProps.amp = obj.amplitude;
            obj.fitProps.bkg = obj.background;
            obj.fitProps.std = obj.stdev;
            obj.fitProps.ampUb = AmpUb;
            obj.fitProps.ampLb = AmpLb;
            obj.fitProps.bkgUb = BkgUb;
            obj.fitProps.bkgLb = BkgLb;
            obj.fitProps.stdUb = SigUb;
            obj.fitProps.stdLb = SigLb;
            obj.fitProps.xCoef = xP;
            obj.fitProps.yCoef = yP;
            obj.fitProps.fixStartPoint = fixStartPoint;
            
            v0 = [ obj.background, obj.amplitude, obj.stdev, xP, yP];
            
            % Create swarm with varying lengths (also useful with
            % lsqnonlin)
            nSwarm = 25;
            matSwarm = repmat(v0, nSwarm, 1);
            tRange = linspace(0.75, 5, nSwarm);
            tOld = linspace(0, 1);

            % get data points from polynomial coefficients
            sIdx = 4;
            eIdx = obj.polyOrder;
            for jSwarm = 1 : nSwarm
                tNew = linspace(0, tRange(jSwarm) );

                xNew = polyval( obj.estimatedCoef(1,:), tNew);
                yNew = polyval( obj.estimatedCoef(2,:), tNew);

                xcf = polyfit( tOld, xNew, obj.polyOrder);
                ycf = polyfit( tOld, yNew, obj.polyOrder);

                if fixStartPoint
                    matSwarm( jSwarm, sIdx : sIdx+eIdx-1) = xcf(1:end-1);
                    matSwarm( jSwarm, sIdx+eIdx: end) = ycf( 1: end-1);
                else
                    matSwarm( jSwarm, sIdx : sIdx+eIdx) = xcf;
                    matSwarm( jSwarm, sIdx+eIdx+1 : end) = ycf; 
                end
            end

            % Create upper and lower bounds for particle swarm
            maxM = max( matSwarm, [], 1);
            minM = min( matSwarm, [], 1);
            matUb = maxM + 0.5*std( matSwarm, 0, 1);
            matLb = minM - 0.5*std( matSwarm, 0, 1);
            if fixStartPoint
                xP_ub = matUb(sIdx : sIdx+eIdx-1);
                yP_ub = matUb(sIdx+eIdx: end);
                xP_lb = matLb(sIdx : sIdx+eIdx-1);
                yP_lb = matLb(sIdx+eIdx: end);
            else
                xP_ub = matUb(sIdx : sIdx+eIdx);
                yP_ub = matUb(sIdx+eIdx+1: end);
                xP_lb = matLb(sIdx : sIdx+eIdx);
                yP_lb = matLb(sIdx+eIdx+1: end);
            end
            obj.fitProps.xCoefUb = xP_ub;
            obj.fitProps.xCoefLb = xP_lb;
            obj.fitProps.yCoefUb = yP_ub;
            obj.fitProps.yCoefLb = yP_lb;
            
            ub = [ BkgUb, AmpUb, SigUb, xP_ub, yP_ub];
            lb = [ BkgLb, AmpLb, SigLb, xP_lb, yP_lb];
            
            obj.fitProps.ub = ub;
            obj.fitProps.lb = lb;
            obj.fitProps.v0 = v0;
            
            FixedParams.PolyOrder = obj.polyOrder;
            FixedParams.StartingPointFixed = obj.fitProps.fixStartPoint;
            FixedParams.XCoefEnd = obj.estimatedCoef(1,end);
            FixedParams.YCoefEnd = obj.estimatedCoef(2,end);
            FixedParams.Mask = 1.0 * logical(obj.source);
            Image2D = obj.source;
            
            if strcmp( fitMethod, 'lsqnonlin')
                obj.fitProps.algorithm = 'lsqnonlin';
                obj.fitProps.opts = optimoptions( @lsqnonlin, 'MaxFunEvals', 1e4, 'OptimalityTolerance', 1e-7, ...
                    'MaxIter', 25, 'TolFun', 1e-6, 'FiniteDifferenceStepSize', 1e-5, ...
                    'FiniteDifferenceType', 'central', 'StepTolerance', 1e-8, 'display', 'iter', 'OutputFcn',@LSQNONLINplotbestf_SA);
                
            elseif strcmp( fitMethod, 'particleswarm')
            
                obj.fitProps.algorithm = 'particleswarm';
                obj.fitProps.opts = optimoptions( @particleswarm, 'ObjectiveLimit', 0, 'SwarmSize',nSwarm,...
                    'InitialSwarmMatrix', matSwarm, 'OutputFcn',@PSplotbestf_SA, ...
                    'MaxIterations', 100, 'FunctionTolerance', 1e-3, 'display', 'iter', ...
                    'UseVectorized', true, 'MinNeighborsFraction',0.25);
                obj.fitProps.swarm = matSwarm;
                
            end
            
            % PSplotbestf_SA {{{
            function stop = PSplotbestf_SA(optimValues,state)
                %PSWPLOTBESTF Plot best function value.
                %
                %   STOP = PSWPLOTBESTF(OPTIMVALUES, STATE) plots OPTIMVALUES.BESTFVAL
                %   against OPTIMVALUES.ITERATION. This function is called from
                %   PARTICLESWARM with the following inputs:
                %
                %   OPTIMVALUES: Information after the current local solver call.
                %          funccount: number of function evaluations
                %              bestx: best solution found so far
                %           bestfval: function value at bestx
                %          iteration: iteration number
                %           meanfval: average function value of swarm particles
                %    stalliterations: number of iterations since improvement in the 
                %                     objective function value stopped
                %              swarm: the position of the swarm particles
                %         swarmfvals: objective function value of swarm particles
                % 
                %   STATE: Current state in which plot function is called. 
                %          Possible values are:
                %             init: initialization state 
                %             iter: iteration state 
                %             done: final state
                %
                %   STOP: A boolean to stop the algorithm.

                % Initialize stop boolean to false.
                stop = false;
                numPixX = size(FixedParams.Mask,1);
                numPixY = size(FixedParams.Mask,2);
                switch state
                    case 'init'

                        t = linspace(0,1); % paramateric variable
                        startIdx = 4;
                        gap = ( length( optimValues.bestx) - startIdx+1)/2;
                        if FixedParams.StartingPointFixed
                            xcurr = polyval( [ optimValues.bestx( startIdx : startIdx+ gap-1), FixedParams.XCoefEnd], t); % current curve coordinates
                            ycurr = polyval( [ optimValues.bestx( startIdx+ gap : end), FixedParams.YCoefEnd], t); % current curve coordinates
                        else
                            xcurr = polyval( optimValues.bestx( startIdx : startIdx+ gap-1), t); % current curve coordinates
                            ycurr = polyval( optimValues.bestx( startIdx+ gap : end), t); % current curve coordinates
                        end

                        figure('Name', ['particleswarm_', num2str(obj.id)], 'NumberTitle', 'off')
                        posOld = get(gcf, 'Position');
                        set(gcf, 'position', [-10000 posOld(2:end)], 'WindowState', 'maximized');
                        subplot(121)
                        plotBest = plot(optimValues.iteration,optimValues.bestfval, '--b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
                        set(plotBest,'Tag','psoplotbestf');
                        xlabel('Iteration','interp','none');
                        ylabel('Function value','interp','none')
                        title(sprintf('Best Function Value: %g',optimValues.bestfval),'interp','none');
                        set(gca, 'FontSize', 14)
                        grid minor; grid on

                        subplot(122)
                        imagesc([Image2D]); axis equal; colormap gray; hold on;
%                         imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow

                    case 'iter'

                        subplot(121)
                        plotBest = findobj(get(gca,'Children'),'Tag','psoplotbestf');
                        newX = [get(plotBest,'Xdata') optimValues.iteration];
                        newY = [get(plotBest,'Ydata') optimValues.bestfval];
                        set(plotBest,'Xdata',newX, 'Ydata',newY);
                        set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimValues.bestfval));
                        grid minor; grid on

                        t = linspace(0,1); % paramateric variable
                        startIdx = 4;
                        gap = ( length( optimValues.bestx) - startIdx+1)/2;
                        if FixedParams.StartingPointFixed
                            xcurr = polyval( [ optimValues.bestx( startIdx : startIdx+ gap-1), FixedParams.XCoefEnd], t); % current curve coordinates
                            ycurr = polyval( [ optimValues.bestx( startIdx+ gap : end), FixedParams.YCoefEnd], t); % current curve coordinates
                        else
                            xcurr = polyval( optimValues.bestx( startIdx : startIdx+ gap-1), t); % current curve coordinates
                            ycurr = polyval( optimValues.bestx( startIdx+ gap : end), t); % current curve coordinates
                        end
                        
                        subplot(122)
                        imagesc([Image2D]); axis equal; colormap gray; hold on;
%                         imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow

                    case 'done'
                        % No clean up tasks required for this plot function.        
                end    


                end
                % }}}

                % LSQNONLINplotbestf_SA {{{
                function stop = LSQNONLINplotbestf_SA(x,optimValues,state)
                %PSWPLOTBESTF Plot best function value.
                %
                %   STOP = PSWPLOTBESTF(OPTIMVALUES, STATE) plots OPTIMVALUES.BESTFVAL
                %   against OPTIMVALUES.ITERATION. This function is called from
                %   PARTICLESWARM with the following inputs:
                %
                %   OPTIMVALUES: Information after the current local solver call.
                %          funccount: number of function evaluations
                %              bestx: best solution found so far
                %           bestfval: function value at bestx
                %          iteration: iteration number
                %           meanfval: average function value of swarm particles
                %    stalliterations: number of iterations since improvement in the 
                %                     objective function value stopped
                %              swarm: the position of the swarm particles
                %         swarmfvals: objective function value of swarm particles
                % 
                %   STATE: Current state in which plot function is called. 
                %          Possible values are:
                %             init: initialization state 
                %             iter: iteration state 
                %             done: final state
                %
                %   STOP: A boolean to stop the algorithm.
                %
                %   See also PARTICLESWARM

                %   Copyright 2014 The MathWorks, Inc.

                % Initialize stop boolean to false.
                stop = false;
                numPixX = size(FixedParams.Mask,1);
                numPixY = size(FixedParams.Mask,2);
                switch state
                    case 'init'

                        t = linspace(0,1); % paramateric variable
                        startIdx = 4;
                        gap = ( length( x) - startIdx+1)/2;
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
%                         imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.5], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow

                    case 'iter'

                        subplot(121)
                        plotBest = findobj(get(gca,'Children'),'Tag','psoplotbestf');
                        newX = [get(plotBest,'Xdata') optimValues.iteration];
                        newY = [get(plotBest,'Ydata') optimValues.resnorm ];
                        set(plotBest,'Xdata',newX, 'Ydata',newY);
                        set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimValues.resnorm) );
                        grid minor; grid on

                        t = linspace(0,1); % paramateric variable
                        startIdx = 4;
                        gap = ( length( x) - startIdx+1)/2;
                        if FixedParams.StartingPointFixed
                            xcurr = polyval( [ x( startIdx : startIdx+ gap-1), FixedParams.XCoefEnd], t); % current curve coordinates
                            ycurr = polyval( [ x( startIdx+ gap : end), FixedParams.YCoefEnd], t); % current curve coordinates
                        else
                            xcurr = polyval( x( startIdx : startIdx+ gap-1), t); % current curve coordinates
                            ycurr = polyval( x( startIdx+ gap : end), t); % current curve coordinates
                        end

                        subplot(122)
                        imagesc([Image2D]); axis equal; colormap gray; hold on;
%                         imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
                        plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
                        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY], 'XTick', [], 'YTick', []);
                        title('Best Curve'); set(gca, 'FontSize', 14)

                        drawnow

                    case 'done'
                        % No clean up tasks required for this plot function.        
                end    


                end
                % }}}
            
        end
        % }}}
        
        % fitMTCurves {{{
        function obj = fitMTCurves( obj)
            
            FixedParams.PolyOrder = obj.polyOrder;
            FixedParams.StartingPointFixed = obj.fitProps.fixStartPoint;
            FixedParams.XCoefEnd = obj.estimatedCoef(1,end);
            FixedParams.YCoefEnd = obj.estimatedCoef(2,end);
            FixedParams.Mask = 1.0 * logical(obj.source);
            Image2D = obj.source;
            
            
            if strcmp( obj.fitProps.algorithm, 'lsqnonlin')
                
                [vfit,resnorm,residual,exitflag,~,~,jacobian] = lsqnonlin( @ErrorFcn, obj.fitProps.v0, obj.fitProps.lb, obj.fitProps.ub, obj.fitProps.opts);
                
            elseif strcmp( obj.fitProps.algorithm, 'particleswarm')
                
                [vfit,resnorm,exitflag,output] = particleswarm( @ErrorFcnPS, length(obj.fitProps.v0), obj.fitProps.lb, obj.fitProps.ub, obj.fitProps.opts);
                
            end
            disp('Mission Accomplished')
            
            
            pStructFin.bkg = vfit(1);
            pStructFin.amp = vfit(2);
            pStructFin.sigma = vfit(3);
            if FixedParams.StartingPointFixed
                pStructFin.XCoef = [ vfit( 4: 4+obj.polyOrder-1), obj.estimatedCoef(1,end)];
                pStructFin.YCoef = [ vfit( 4+obj.polyOrder :end), obj.estimatedCoef(2,end)];
            else
                pStructFin.XCoef = vfit( 4: 4+obj.polyOrder-1);
                pStructFin.YCoef = vfit( 4+ obj.polyOrder :end);
            end
            
            obj.fitProps.vfit = vfit;
            obj.fitProps.structFinal = pStructFin;
            
            if pStructFin.amp < pStructFin.bkg
                obj.dead = 1;
            end
            
            
                function err = ErrorFcn( p)

                   err = GaussianCurveMaker2D( p, Image2D, FixedParams) - Image2D;
            %        err = sum ( err(:).^2 );
                   err = err(:);

                end

                function err = ErrorFcnPS( p)

                   errMat = GaussianCurveMaker2D_PSvectorized( p, Image2D, FixedParams) - repmat( Image2D, 1, 1, nSwarm);
                    err = squeeze( sum( sum(errMat.^2, 1), 2) );
                end

                % GaussianCurveMaker2D {{{
                function imPlane = GaussianCurveMaker2D( p, Image2D, FixedParams)

                    order = FixedParams.PolyOrder;
                    Bkg = p(1); Amp = p(2); Sig = p(3);

                    if FixedParams.StartingPointFixed == 1
                        XCoef = [ p(4: 4+order-1), FixedParams.XCoefEnd ];
                        YCoef = [ p(4+order: end), FixedParams.YCoefEnd ];
                    elseif FixedParams.StartingPointFixed == 0
                        XCoef = p(4: 4+order-1);
                        YCoef = p(4+order: end);
                    end

                    % Initialize the image volume
                    imPlane = 0*Image2D;
                    numPixX = size(imPlane, 1);
                    numPixY = size(imPlane, 2);

                    % load the meshgrids
                    [ yGrid, xGrid ] = meshgrid( 1 : numPixY , 1 : numPixX);
                    xV = xGrid(:); yV = yGrid(:);
                    idx = 1 : numel( imPlane);

                    % create parametric coord
                    Tvec = linspace(0,1);
                    Tvec( end) = [];

                    % Now for each value of the parameter, we'll do a gauss
                    % quadrature numerical integration 
                    DeltaT = median( diff(Tvec) );

                    % get the offsets for gauss quadrature
                    poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
                    poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;

                    % for speed.
                    Tsort = sort([  Tvec + poff1,  Tvec + poff2] );
                    xLoc = polyval( XCoef, Tsort )';
                    yLoc = polyval( YCoef, Tsort )';

                    dParam = 1;
                    for jk = 1 : length( xLoc) -1
                        dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2)];
                    end

                    fitX = fit( dParam', xLoc, 'linearinterp');
                    fitY = fit( dParam', yLoc, 'linearinterp');
                    dParamNew = linspace( dParam(1), dParam(end), length(dParam) );
                    xDat = fitX( dParamNew);
                    yDat = fitY( dParamNew);

                    PtGaussLoc = [ xDat , yDat ];


                    % Compute the approximate integral over the spline points by
                    % multiplying gaussians in x, y and z and adding for each
                    % different value of the parameter.
                    Conv2 = NumericalConv2( Sig, PtGaussLoc, xV, yV, imPlane, idx );
                    imFeatures = ( Amp - Bkg) * mat2gray(Conv2);
                    % set very small values to 0
                    imFeatures ( imFeatures < 0.1*max( imFeatures(:) ) ) = 0;

                    imBkg = (imPlane + Bkg) .* FixedParams.Mask;
                    imPlane = imBkg + ( Amp - Bkg) * mat2gray( imFeatures);

                    % Penalize lsqnonlin if it draws a curve outside the mask region
                    imCheck = imPlane .* imcomplement( FixedParams.Mask);
                    if any( imCheck(:) > 0.9*max(imFeatures(:) ) )
                        imPlane = imPlane + Inf.*imcomplement( FixedParams.Mask);
                    end

                    function imConv2 = NumericalConv2( Sig, PtGaussLoc, xV, yV, imPlane, idx )

                        % Create copies of image volume
                        imConv2 = imPlane .* 0;

                        NumGauss = size( PtGaussLoc, 1);

                        CoordsX = yV .* ones( 1, NumGauss);
                        CoordsY = xV .* ones( 1, NumGauss);

                        xx = PtGaussLoc( :, 1);
                        yy = PtGaussLoc( :, 2);

                        s1 = Sig;
                        s2 = Sig;

                        Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
                        Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);

                        imConv2( idx') = sum( Xvals .* Yvals , 2);

                    end


                end
                % }}}

                % GaussianCurveMaker2D_PSvectorized {{{
                function imPlaneVec = GaussianCurveMaker2D_PSvectorized( p, Image2D, FixedParams)

                    order = FixedParams.PolyOrder;
                    Bkg = p(:,1); Amp = p(:,2); Sig = p(:,3);
                    numSwarm = size(p,1);

                    if FixedParams.StartingPointFixed == 1
                        XCoef = [ p(:, 4: 4+order-1), repmat(FixedParams.XCoefEnd, size(p,1), 1) ];
                        YCoef = [ p(:, 4+order: end), repmat(FixedParams.YCoefEnd, size(p,1), 1) ];
                    elseif FixedParams.StartingPointFixed == 0
                        XCoef = p(:, 4: 4+order-1);
                        YCoef = p(:, 4+order: end);
                    end

                    % Initialize the image volume
                    imPlane = 0*Image2D;
                    numPixX = size(imPlane, 1);
                    numPixY = size(imPlane, 2);

                    % load the meshgrids
                    [ yGrid, xGrid ] = meshgrid( 1 : numPixY , 1 : numPixX);
                    xV = xGrid(:); yV = yGrid(:);
                    idx = 1 : numel( imPlane);

                    % create parametric coord
                    Tvec = linspace(0,1);
                    Tvec( end) = [];

                    % Now for each value of the parameter, we'll do a gauss
                    % quadrature numerical integration 
                    DeltaT = median( diff(Tvec) );

                    % get the offsets for gauss quadrature
                    poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
                    poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;

                    % for speed.
                    Tsort = sort([  Tvec + poff1,  Tvec + poff2] );

                    xLoc = zeros(length(Tsort), numSwarm); yLoc = xLoc;
                    for jSwarm = 1 : numSwarm
                        xLoc(:,jSwarm) = polyval( XCoef(jSwarm,:), Tsort )';
                        yLoc(:,jSwarm) = polyval( YCoef(jSwarm,:), Tsort )';
                    end

                    dParam = ones(numSwarm, 1);
                    for jk = 1 : size( xLoc,1) -1
                        dParam = [ dParam(:,:), dParam(:,end) + ( sqrt( diff( xLoc( jk:jk+1, :), 1, 1).^2 + diff( yLoc( jk:jk+1, : ), 1, 1 ).^2) )'];
                    end

                    xDat = 0*xLoc; yDat = 0*yLoc;
                    fitXvec = cell( numSwarm, 1); fitYvec = fitXvec;
                    dParamNew = zeros(size(dParam'));
                    PtGaussLoc = zeros( size( xLoc, 1), 2, numSwarm);
                    for jSwarm = 1 : numSwarm
                        fitXvec{jSwarm} = fit( dParam(jSwarm, :)', xLoc(:,jSwarm), 'linearinterp');
                        fitYvec{jSwarm} = fit( dParam(jSwarm, :)', yLoc(:,jSwarm), 'linearinterp');

                        dParamNew(:,jSwarm) = linspace( dParam(jSwarm, 1), dParam(jSwarm,end), length( dParam(jSwarm,:) ) );

                        xDat(:,jSwarm) = fitXvec{jSwarm}( dParamNew( :,jSwarm) );
                        yDat(:,jSwarm) = fitYvec{jSwarm}( dParamNew( :,jSwarm) );

                        PtGaussLoc(:,:,jSwarm) = [ xDat(:,jSwarm) , yDat(:,jSwarm) ];

                    end

                    % Compute the approximate integral over the spline points by
                    % multiplying gaussians in x, y and z and adding for each
                    % different value of the parameter.
                    Conv2Vec = NumericalConv2_PSvectorized( Sig, PtGaussLoc, xV, yV, imPlane, idx );
                    imPlaneVec = 0*Conv2Vec;
                    for jSwarm = 1 : numSwarm
                        imPlane = zeros( size( imPlaneVec, 1), size( imPlaneVec, 2) );
                        imFeatures = ( Amp(jSwarm) - Bkg(jSwarm)) * mat2gray( Conv2Vec(:,:,jSwarm) );
                        % set very small values to 0
                        imFeatures ( imFeatures < 0.1*max( imFeatures(:) ) ) = 0;

                        imBkg = (imPlane + Bkg(jSwarm) ) .* FixedParams.Mask;

                        imPlane = imBkg + ( Amp(jSwarm) - Bkg(jSwarm)) * mat2gray( imFeatures);

                        % Penalize lsqnonlin if it draws a curve outside the mask region
                        imCheck = imPlane .* imcomplement( FixedParams.Mask);
                        if any( imCheck(:) > 0.9*max(imFeatures(:) ) )
                            imPlane = imPlane + Inf.*imcomplement( FixedParams.Mask);
                            imPlane( isnan(imPlane(:) ) ) = 0;
                        end
                        imPlaneVec(:,:,jSwarm) = imPlane;

                        if any(isnan(imPlaneVec(:) ) )
                            stopHere = 1;
                        end
                    end

                        function imConvVec = NumericalConv2_PSvectorized( Sig, PtGaussLoc, xV, yV, imPlane, idx )

                            NumGauss = size( PtGaussLoc, 1);
                            NumSwarm = size( PtGaussLoc, 3);

                            % Create copies of image volume
                            imConv2 = imPlane .* 0;
                            imConvVec = repmat(imConv2, 1, 1, NumSwarm);
                            for jSwarm = 1 : numSwarm

                                CoordsX = yV .* ones( 1, NumGauss, 1);
                                CoordsY = xV .* ones( 1, NumGauss, 1);

                                xx = PtGaussLoc( :, 1, jSwarm);
                                yy = PtGaussLoc( :, 2, jSwarm);

                                s1 = Sig(jSwarm);
                                s2 = Sig(jSwarm);

                                Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
                                Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);

                                imConv2( idx') = sum( Xvals .* Yvals , 2);

                                imConvVec(:,:,jSwarm) = imConv2;

                            end

                        end


                end
                
                % }}}

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

