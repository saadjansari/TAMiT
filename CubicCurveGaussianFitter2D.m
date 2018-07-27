function pStructFin = CubicCurveGaussianFitter2D( pStruct, Image2D, algorithm)

   
    
    %% Load and Initialize values
           
    % Fix starting point
    StartingPointFixed = 1;
    
    XCoefInit = pStruct.XCoef;
    YCoefInit = pStruct.YCoef;
    AmpInit = pStruct.amp;
    SigInit = pStruct.sigma;
    BkgInit = pStruct.bkg;
    
    numPixX = size(Image2D, 1);
    numPixY = size(Image2D, 2);
    
    if StartingPointFixed
        FixedParams.XCoefEnd = XCoefInit(end);
        FixedParams.YCoefEnd = YCoefInit(end);
        XCoefInit( end) = [];
        YCoefInit( end) = [];
    end
    FixedParams.StartingPointFixed = StartingPointFixed;
    FixedParams.PolyOrder = length( pStruct.XCoef) - 1;
    FixedParams.Mask = logical( Image2D);
    
    
    
    %% Set upper and lower bounds for parameters.
    
    AmpD = 0.01;
    SigD = 0;
    BkgD = 0.01;
    AmpUb = AmpInit + AmpD; AmpLb = AmpInit - AmpD;
    SigUb = SigInit + SigD; SigLb = SigInit - SigD; 
    BkgUb = BkgInit + BkgD; BkgLb = BkgInit - BkgD;
    
    % create upper lower bound vector
    XCoefUb = XCoefInit; XCoefUb( XCoefUb < 0) = 0; XCoefUb( XCoefUb > 0) = +Inf;
    XCoefLb = XCoefInit; XCoefLb( XCoefLb < 0) = -Inf; XCoefLb( XCoefLb > 0) = 0;
    YCoefUb = YCoefInit; YCoefUb( YCoefUb < 0) = 0; YCoefUb( YCoefUb > 0) = +Inf;
    YCoefLb = YCoefInit; YCoefLb( YCoefLb < 0) = -Inf; YCoefLb( YCoefLb > 0) = 0;
    if StartingPointFixed == 0
       
        XCoefUb(end) = numPixX; YCoefUb(end) = numPixY;
        XCoefLb(end) = 1; YCoefLb(end) = 1;
        
    end
    v0 = [ BkgInit, AmpInit, SigInit, XCoefInit, YCoefInit];
    ub = [ BkgUb, AmpUb, SigUb, XCoefUb, YCoefUb];
    lb = [ BkgLb, AmpLb, SigLb, XCoefLb, YCoefLb];
    
%     coefUb = ones(1, length(XCoefInit) ) * Inf;
%     coefLb = ones(1, length(XCoefInit) ) * -Inf;
%     
    % Finally construct upper and lower bounds
%     v0 = [ BkgInit, AmpInit, SigInit, XCoefInit, YCoefInit];
%     ub = [ BkgUb, AmpUb, SigUb, coefUb, coefUb];
%     lb = [ BkgLb, AmpLb, SigLb, coefLb, coefLb];
    
    % Create particle swarm
    nSwarm = 25;
    matSwarm = repmat(v0, nSwarm, 1);
    tRange = linspace(0.75, 1.25, nSwarm);
    tOld = linspace(0, 1);
    
%     h2 = figure; 
    % get data points from polynomial coefficients
    for jSwarm = 1 : nSwarm
        tNew = linspace(0, tRange(jSwarm) );

        xNew = polyval( pStruct.XCoef, tNew);
        yNew = polyval( pStruct.YCoef, tNew);

        xcf = polyfit( tOld, xNew, FixedParams.PolyOrder);
        ycf = polyfit( tOld, yNew, FixedParams.PolyOrder);
        
%         figure(h2); imagesc(Image2D); hold on; plot( polyval( xcf, tOld), polyval( ycf, tOld), 'r-', 'LineWidth', 5); title(sprintf('t = %.1f',tRange(jSwarm)))
%         pause(1)
        
        sIdx = 4;
        eIdx = FixedParams.PolyOrder;
        if StartingPointFixed
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
    ubPS = [ BkgUb, AmpUb, SigUb, matUb(4: end)];
    lbPS = [ BkgLb, AmpLb, SigLb, matLb(4: end)];
    
    if strcmp(algorithm, 'ps')
        optsPS = optimoptions( @particleswarm, 'ObjectiveLimit', 0, 'SwarmSize',nSwarm,...
            'InitialSwarmMatrix', matSwarm, 'OutputFcn',@PSplotbestf_SA, ...
            'MaxIterations', 100, 'FunctionTolerance', 1e-3, 'display', 'iter', ...
            'UseVectorized', true, 'MinNeighborsFraction',0.25);
        vfit = particleswarm( @ErrorFcnPS, length(v0), lbPS, ubPS, optsPS);
    elseif strcmp(algorithm, 'lsqnonlin')
        % define fitting options
        opts = optimoptions( @lsqnonlin, 'MaxFunEvals', 1e4, 'OptimalityTolerance', 1e-7, ...
            'MaxIter', 100, 'TolFun', 1e-6, 'FiniteDifferenceStepSize', 1e-5, ...
            'FiniteDifferenceType', 'central', 'StepTolerance', 1e-8, 'display', 'iter', 'OutputFcn',@LSQNONLINplotbestf_SA);

        % Run lsqnonlin
        [vfit,resnorm,residual,exitflag,~,~,jacobian] = lsqnonlin( @ErrorFcn, v0, lbPS, ubPS, opts);
    end
%     disp(['ExitFlag = ', num2str(exitflag)])
    disp('Mission Accomplished')
    
    numTubes = 2;
    colors = distinguishable_colors( numTubes, {'w', 'k'} );
    
    pStructFin.bkg = vfit(1);
    pStructFin.amp = vfit(2);
    pStructFin.sigma = vfit(3);
    if FixedParams.StartingPointFixed
        pStructFin.XCoef = [ vfit( 4: 4+FixedParams.PolyOrder-1), FixedParams.XCoefEnd];
        pStructFin.YCoef = [ vfit( 4+FixedParams.PolyOrder :end), FixedParams.YCoefEnd];
    else
        pStructFin.XCoef = vfit( 4: 4+FixedParams.PolyOrder-1);
        pStructFin.YCoef = vfit( 4+ FixedParams.PolyOrder :end);
    end
    
    compareFitPlot = 0;
    % Compare fitting
    if compareFitPlot
        
        t = linspace(0,1); % paramateric variable
        xinit = polyval( pStruct.XCoef, t); % init and final curve coordinates
        yinit = polyval( pStruct.YCoef, t);
        xfin = polyval( pStructFin.XCoef, t);
        yfin = polyval( pStructFin.YCoef, t);

        % plot vars
        trans = 0.5;
        colInit = [ colors(1, :), trans];
        colFin = [ colors(2, :), trans];
        minVal = min( Image2D(:) );
        maxVal = max( Image2D(:) );
        clim = [minVal, maxVal];
        lw = 5;
        
        figure;
        subplot(131); imagesc( Image2D, clim); colormap gray; axis equal; hold on; colorbar;
        plot( xinit, yinit, 'color', colFin);
        plot( xfin, yfin, 'color', colFin); hold off; title('Fitted Curve')
        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY]);        
        subplot(132); imagesc( GaussianCurveMaker2D( v0, Image2D, FixedParams), clim ); colormap gray; axis equal; colorbar
        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY]); title('Estimated Intensity')
        subplot(133); imagesc( GaussianCurveMaker2D( vfit, Image2D, FixedParams), clim ); colormap gray; axis equal; colorbar;
        set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY]); title('Fitted Intensity'); 
%         set(gcf, 'pos', get(0, 'ScreenSize') );
        
        figure; imagesc( Image2D, clim); colormap gray; axis equal; hold on; 
        plot( xinit, yinit, 'Color', colInit, 'LineWidth', lw); 
        plot( xfin, yfin, 'Color', colFin, 'LineWidth', lw); 
        % create legend entry
        LH(1) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colInit); L{1} = 'Estimated Line';
        LH(2) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colFin); L{2} = 'Fitted Line';
        legend(LH, L); hold off; set(gca, 'FontSize', 15, 'xlim', [1 numPixX], 'ylim', [1 numPixY]);
        
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
        
%         xLoc = sort( [V1X ; V2X]);
%         yLoc = sort( [V1Y ; V2Y]);
%         PtGaussLoc = [ xLoc , yLoc ];

%         PtGaussLoc = [ xLoc , yLoc ];
        
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
%         imPlaneRaw = ( Amp - Bkg) * mat2gray(Conv2) + Bkg;
%         imPlane = imPlaneRaw .* FixedParams.Mask;
        
        
%         imPenalize = ( imcomplement( FixedParams.Mask) .* Conv2);
%         if any( imPenalize(:) > 0.9*max( Conv2(:) ) )
%             imPlane(:) = Inf;
%         end
%         figure; subplot(131); imagesc(Image2D)
%         subplot(132); imagesc( imPlaneRaw);
%         subplot(133); imagesc( imPenalize);

%         figure; 
%         subplot(121); imagesc(Image2D); colormap gray; axis equal; hold on; plot(xLoc, yLoc, 'r-', 'LineWidth', 3); hold off
%         subplot(122); imagesc(imPlane); colormap gray; axis equal; 
        
        
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
        
%         % Penalize lsqnonlin if it draws a curve outside the mask region
%         imCheck = imPlane .* imcomplement( FixedParams.Mask);
%         if any( imCheck(:) > 0.9*max(imFeatures(:) ) )
%             imPlane = imPlane + Inf.*imcomplement( FixedParams.Mask);
%         end
        
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
            
%                 tic
%                 % break it off into 4 pieces
%                 numPiece = 4;
%                 len = length(idx);
%                 off = mod( len, numPiece );
%                 pieceSize = (len-off)/numPiece;
%                 
%                 for jChunk = 1: numPiece
%                     
%                     indices = ( (jChunk-1)*pieceSize +1 : jChunk*pieceSize + off*(jChunk==numPiece) );
%                     idxChunk = idx(  indices );
%                     
%                     CoordsX = yV( indices) .* ones( 1, NumGauss, NumSwarm);
%                     CoordsY = xV( indices) .* ones( 1, NumGauss, NumSwarm);
% 
%                     xx = permute( PtGaussLoc( :, 1, :), [2 1 3]);
%                     yy = permute( PtGaussLoc( :, 2, :), [2 1 3]);
%                     
%                     s1 = zeros(1,1,NumSwarm); s2 = s1;
%                     s1(1,1,:) = Sig;
%                     s2(1,1,:) = Sig;
% 
%                     Xvals = exp( -( xx - CoordsX ).^2 ./ (2 * s1.^2) ) ./ ( 2 * s1.^2 * pi).^(1/2);
%                     Yvals = exp( -( yy - CoordsY ).^2 ./ (2 * s2.^2) ) ./ ( 2 * s2.^2 * pi).^(1/2);
% 
%                     for jSwarm = 1 : NumSwarm
%                         imConv2( idxChunk') = sum( Xvals(:,:,jSwarm) .* Yvals(:,:,jSwarm) , 2);
% 
%                         imConvVec(:,:,jSwarm) = imConv2;
%                     end
%                     
%                     
%                 end
%                 toc
                
            end

        
    end

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
    %
    %   See also PARTICLESWARM

    %   Copyright 2014 The MathWorks, Inc.

    % Initialize stop boolean to false.
    stop = false;
    switch state
        case 'init'
%             plotBest = plot(optimValues.iteration,optimValues.bestfval, '-b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
%             set(plotBest,'Tag','psoplotbestf');
%             xlabel('Iteration','interp','none');
%             ylabel('Function value','interp','none')
%             title(sprintf('Best Function Value: %g',optimValues.bestfval),'interp','none');
%             set(gca, 'FontSize', 14)
%             grid on
            
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
%             figure('Name', 'particle swarm best curve')
%             imagesc(Image2D); axis equal; colormap gray; hold on;
%             plot( xcurr, ycurr, 'Color', 'r', 'LineWidth', 4); 
%             set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY]);
            
            figure('Name', 'particle swarm', 'NumberTitle', 'off')
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
            imagesc([Image2D; GaussianCurveMaker2D( optimValues.bestx, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
            plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
            set(gca, 'xlim', [1 numPixX], 'ylim', [1 2*numPixY], 'XTick', [], 'YTick', []);
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
            imagesc([Image2D; GaussianCurveMaker2D( optimValues.bestx, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
            plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
            set(gca, 'xlim', [1 numPixX], 'ylim', [1 2*numPixY], 'XTick', [], 'YTick', []);
            title('Best Curve'); set(gca, 'FontSize', 14)
            
            drawnow
            
        case 'done'
            % No clean up tasks required for this plot function.        
    end    


    end



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
    switch state
        case 'init'
%             plotBest = plot(optimValues.iteration,optimValues.bestfval, '-b', 'Marker', '*', 'LineWidth', 3, 'MarkerSize', 10);
%             set(plotBest,'Tag','psoplotbestf');
%             xlabel('Iteration','interp','none');
%             ylabel('Function value','interp','none')
%             title(sprintf('Best Function Value: %g',optimValues.bestfval),'interp','none');
%             set(gca, 'FontSize', 14)
%             grid on
            
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
%             figure('Name', 'particle swarm best curve')
%             imagesc(Image2D); axis equal; colormap gray; hold on;
%             plot( xcurr, ycurr, 'Color', 'r', 'LineWidth', 4); 
%             set(gca, 'xlim', [1 numPixX], 'ylim', [1 numPixY]);
            
            figure('Name', 'LSQNONLIN', 'NumberTitle', 'off');
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
            imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
            plot( xcurr, ycurr, 'Color', [1 0 0 0.5], 'LineWidth', 6); hold off
            set(gca, 'xlim', [1 numPixX], 'ylim', [1 2*numPixY], 'XTick', [], 'YTick', []);
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
            imagesc([Image2D; GaussianCurveMaker2D( x, Image2D, FixedParams) ]); axis equal; colormap gray; hold on;
            plot( xcurr, ycurr, 'Color', [1 0 0 0.7], 'LineWidth', 6); hold off
            set(gca, 'xlim', [1 numPixX], 'ylim', [1 2*numPixY], 'XTick', [], 'YTick', []);
            title('Best Curve'); set(gca, 'FontSize', 14)
            
            drawnow
            
        case 'done'
            % No clean up tasks required for this plot function.        
    end    


    end


end

