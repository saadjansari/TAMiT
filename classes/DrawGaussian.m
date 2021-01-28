function [imageFeat, error_code, error_amount] = DrawGaussian( sigma, imageIn, ftype, varargin)
% Draw gaussian feature of type ftype
% ftype: 'Spot2', 'Spot3', 'Line2', 'Line3', 'Curve2', 'Curve3'

    opts = parseArgs( imageIn, varargin{:});

    switch ftype
        case 'Spot2'
            [imageFeat, error_code, error_amount] = drawGaussianPoint2D( opts.Pos, sigma, imageIn, opts.Idx, opts.X, opts.Y);
        case 'Spot3'
            [imageFeat, error_code, error_amount] = drawGaussianPoint3D( opts.Pos, sigma, imageIn, opts.Idx, opts.X, opts.Y, opts.Z);
        case 'Line2'
            [imageFeat, error_code, error_amount] = drawGaussianLine2D( opts.PosStart, opts.PosEnd, sigma, imageIn, opts.Idx, opts.X, opts.Y);
        case 'Line3'
            [imageFeat, error_code, error_amount] = drawGaussianLine3D( opts.PosStart, opts.PosEnd, sigma, imageIn, opts.Idx, opts.X, opts.Y, opts.Z);
        case 'Curve2'
            [imageFeat, error_code, error_amount] = drawGaussianCurve2D( opts.Coeff, sigma, imageIn, opts.T, opts.Idx, opts.X, opts.Y);
        case 'Curve3'
            [imageFeat, error_code, error_amount] = drawGaussianCurve3D( opts.Coeff, sigma, imageIn, opts.T, opts.Idx, opts.X, opts.Y, opts.Z);
        case 'Curve3Coords'
            [imageFeat, error_code, error_amount] = drawGaussianCurve3DCoords( opts.Coord, sigma, imageIn, opts.Idx, opts.X, opts.Y, opts.Z);
        case 'Curve2Coords'
            [imageFeat, error_code, error_amount] = drawGaussianCurve2DCoords( opts.Coord, sigma, imageIn, opts.Idx, opts.X, opts.Y);
        otherwise
            error('DrawGaussian: unknown feature type')
    end

    if opts.Plot
        figure;
        if opts.Dim == 2
            imshow( imageFeat);
        elseif opts.Dim == 3
            imshow3( imageFeat);
        end
    end


    % parseArgs {{{
    function opts = parseArgs( imageIn, varargin)

        % Dimension
        defDim = length( size(imageIn));

        % Position for Spot
        defPos = [];

        % Start End Position for Line
        defPosStart = [];
        defPosEnd = [];

        % Coefficients
        defCoeff = [];

        % Coefficients
        defCoord = [];
        
        % Index/XYZ
        defIdx = 1 : numel( imageIn);
        if defDim == 2
            [defy, defx] = ind2sub( size( imageIn), defIdx);
            defx=defx'; defy=defy'; defz = 0*defx;
        elseif defDim == 3
            [defy, defx, defz] = ind2sub( size( imageIn), defIdx);
            defx=defx'; defy=defy';defz=defz';
        end

        % Parameter T for parametric curves
        defT = [0 1];

        % Plot
        defPlot = 0;

        % Input Parser
        p = inputParser;
        addParameter( p, 'Dim', defDim);
        addParameter( p, 'Pos', defPos);
        addParameter( p, 'PosStart', defPosStart);
        addParameter( p, 'PosEnd', defPosEnd);
        addParameter( p, 'Coeff', defCoeff);
        addParameter( p, 'X', defx);
        addParameter( p, 'Y', defy);
        addParameter( p, 'Z', defz);
        addParameter( p, 'Idx', defIdx);
        addParameter( p, 'T', defT);
        addParameter( p, 'Plot', defPlot);
        addParameter( p, 'Coord', defCoord);

        parse( p, varargin{:});
        opts = p.Results;


    end
    % }}}

    % drawGaussianPoint2D {{{
    function [imageSpot, errorCode, error_amount] = drawGaussianPoint2D( pos, sigma, imageIn, idx, x, y)
        % Draws a gaussian 2D point using an analytical framework

        if isempty( pos)
            error('drawGaussianPoint2D: position not provided')
        end

        % Check that everything is in 2D
        if numel( size(imageIn) ) ~= 2 || length(pos)~=2 || length(sigma)~=2
            error('drawGaussianPoint2D: all data must be 2-dimensional')
        end
        
        % Error-checking
        errorCode = 0;
        error_amount = 0;
        if pos(1) > size(imageIn,2) || pos(2) > size(imageIn,1) ...
                || pos(1) < 1 || pos(2) < 1
            error_amount = 1;
            imageSpot = 0*imageIn;
            return
        end

        % Amplitudes:
        ExpX = exp( -0.5*( (x- pos(1))./sigma(1)).^2 );
        ExpY = exp( -0.5*( (y- pos(2))./sigma(2)).^2 );
        IntValues = ExpX .* ExpY;
%         IntValues( isnan( IntValues) ) = min( IntValues(:) );
%         IntValues( IntValues == Inf) = min( IntValues(:) );

        % Initiliaze the volume and set the appropriate indices to these values
        imageSpot = 0*imageIn;
        imageSpot( idx) = IntValues;

    end
    % }}}

    % drawGaussianPoint3D {{{
    function [imageSpot, errorCode, error_amount] = drawGaussianPoint3D( pos, sigma, imageIn, idx, x, y, z)
        % Draws a gaussian point using an analytical framework

        if isempty( pos)
            error('drawGaussianPoint3D: position not provided')
        end
        % Check that everything is in 3D
        if numel( size(imageIn) ) ~= 3 || length(pos)~=3 || length(sigma)~=3
            error('drawGaussianPoint3D: all data must be 3-dimensional')
        end
        % Error-checking
        errorCode = 0;
        error_amount = 0;
        if pos(1) > size(imageIn,2) || pos(2) > size(imageIn,1) ...
                || pos(1) < 1 || pos(2) < 1 
            error_amount = 1;
            imageSpot = 0*imageIn;
            return
        end
        if pos(3) < 1
            error_amount = abs( 1-pos(3));
            errorCode = 1;
        elseif pos(3) > size(imageIn,3)
            error_amount = abs( pos(3) - size(imageIn,3));
            errorCode = size(imageIn,3);
        end

        % Amplitudes:
        ExpX = exp( -0.5*( (x-pos(1))./sigma(1)).^2 );
        ExpY = exp( -0.5*( (y-pos(2))./sigma(2)).^2 );
        ExpZ = exp( -0.5*( (z-pos(3))./sigma(3)).^2 );
        IntValues = ExpX .* ExpY .* ExpZ;
        %IntValues( isnan( IntValues) ) = min( IntValues(:) );
        %IntValues( IntValues == Inf) = min( IntValues(:) );

        % Initialize the volume and set the appropriate indices to these values
        imageSpot = 0*imageIn;
        imageSpot( idx) = IntValues;

    end
    % }}}

    % drawGaussianLine2D {{{
    function [imageLine, errorCode, error_amount] = drawGaussianLine2D( startPos, endPos, sigma, imageIn, idx, x, y)
        % Draws a gaussian straight line using an analytical framework

        if isempty( startPos) || isempty(endPos)
            error('drawGaussianLine2D: start and/or end position not provided')
        end
        % Check that everything is in 2D
        if numel( size(imageIn) ) ~= 2 || length(startPos)~=2 || length(endPos)~=2 || length(sigma)~=2
            error('drawGaussianLine2D: all data must be 2-dimensional')
        end
        errorCode = 0; error_amount = 0;
        
        x0 = startPos(1); y0 = startPos(2);
        x1 = endPos(1); y1 = endPos(2);
        sx = sigma(1); sy = sigma(2);

        % First lets parameterize this line with a parameter t. We'll find the
        % slopes of the line in each spatial dimension along with any offset
        % At t = 0, x = x0, y = y0, z = z0
        % At t = 1, x = x1, y = y1, z = z1
        t0 = 0;
        t1 = 1;

        % find slopes
        mx = (x1 - x0) / ( t1 - t0);
        my = (y1 - y0) / ( t1 - t0);

        % find offsets, these are just the start points at t = 0;
        cx = x0;
        cy = y0;

        % We will integrate a 2D gaussian with respect to the parameter t going
        % from 0 to 1.

        % We've done the analytical integration in mathematica and will use the
        % result here:

        
        ExpVal = - ( (y-cy)*mx + (-x+cx)*my ).^2 ./ (my^2 * sx^2 + mx^2 *sy^2);
        
        Amp = sqrt(pi) * sx * sy ./ (2*sqrt(my^2 *sx^2 + mx^2 * sy^2) );
        
        Erf1 = - erf( (sx^2 *my * (-y + cy) + sy^2 *mx * (-x+cx) ) ./ (sx*sy*sqrt(my^2 *sx^2 + mx^2 * sy^2) ) );
        
        Erf2 = + erf( (sx^2 *my * (-y + cy + my) + sy^2 *mx * (-x+cx + mx) ) ./ (sx*sy*sqrt(my^2 *sx^2 + mx^2 * sy^2) ) );
        
        IntValues = Amp .* exp( ExpVal) .* ( Erf1 + Erf2 );

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

    % drawGaussianLine3D {{{
    function [imageLine, errorCode, error_amount] = drawGaussianLine3D( startPos, endPos, sigma, imageIn, idx, x, y, z)
        % Draws a gaussian straight line using an analytical framework

        if isempty( startPos) || isempty(endPos)
            error('drawGaussianLine2D: start and/or end position not provided')
        end
        % Check that everything is in 3D
        if numel( size(imageIn) ) ~= 3 || length(startPos)~=3 || length(endPos)~=3 || length(sigma)~=3
            error('drawGaussianLine3D: all data must be 3-dimensional')
        end
        
        x0 = startPos(1); y0 = startPos(2); z0 = startPos(3);
        x1 = endPos(1); y1 = endPos(2); z1 = endPos(3);
        sx = sigma(1); sy = sigma(2); sz = sigma(3);
        
        if z0 < 1 || z0 > size(imageIn,3)
            error('3D line start point should not be outside of z region')
        end
        errorCode = 0;
        error_amount = 0;
        
        if z1 < 1
            errorCode = 1;
        elseif z1 > size(imageIn,3)
            errorCode = size(imageIn,3);
        end
        if z1 < 1 || z1 > size(imageIn,3)
            [x1,y1,z1,error_amount] = TrimAndGetErrorZ( startPos,endPos,imageIn);
        end

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
        
        % TrimAndGetErrorZ {{{
        function [X,Y,Z,error_amount] = TrimAndGetErrorZ( startPos, endPos,imageIn)
            
            fun = @(idx) linspace(startPos(idx),endPos(idx),100);
            xx = fun(1); yy = fun(2); zz = fun(3);
            
            if endPos(3) < 1
                idxGood = find( zz < 1, 1, 'first') -1;
            elseif endPos(3) > size(imageIn,3)
                idxGood = find( zz > size(imageIn,3), 1, 'first') -1;
            end
            X = xx( idxGood); Y = yy( idxGood); Z = zz( idxGood);
            error_amount = (100-idxGood)/100;
        end
        % }}}

    end
    % }}}

    % drawGaussianCurve2D {{{
    function imageCurve = drawGaussianCurve2D( coeffs, sigma, imageIn, idx, x,y)
        % Draw a gaussian curve in 3D using numerical integration

        if isempty( coeffs) 
            error('drawGaussianCurve2D: coeffs not provided')
        end
        if numel( size(imageIn) ) ~= 2 || length( coeffs) ~= 2 || length(sigma)~= 2
            error('drawGaussianLine2D: all data must be 2-dimensional')
        end

        if nargin < 4
            idx = 1 : numel( imageIn);
        end
        if nargin < 6
            [y, x] = ind2sub( size( imageIn), idx);
            x=x'; y=y';
        end

        dim = length( size(imageIn) );
        % Create coordinates for numerical integration
        % create parametric coord
        speedUp = 2.0;
        Tvec = linspace(0,1, round(100/speedUp) );
        Tvec( end) = [];

        % Now for each value of the parameter, we'll do a gauss
        % quadrature numerical integration 
        DeltaT = median( diff(Tvec) );

        % get the offsets for gauss quadrature
        poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
        poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;
        Tsort = sort([  Tvec + poff1,  Tvec + poff2] );

        % Coefficients
        XCoef = coeffs{1};
        YCoef = coeffs{2};
        xLoc = polyval( XCoef, Tsort )';
        yLoc = polyval( YCoef, Tsort )';

        % re-parametrize by length
        dParam = 1;
        for jk = 1 : length( xLoc) -1
            dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2)];
        end
        dParamNew = linspace( dParam(1), dParam(end), length(dParam) );

        % fit once again, and find now parametric coordinates
        fitX = fit( dParam', xLoc, 'linearinterp');
        fitY = fit( dParam', yLoc, 'linearinterp');
        xDat = fitX( dParamNew);
        yDat = fitY( dParamNew);

        PtGaussLoc = [ xDat , yDat ];
        imageCurve = NumericalConv( sigma, PtGaussLoc, {x , y}, imageIn, idx );

    end
    % }}}

    % drawGaussianCurve3D {{{
    function [imageCurve, error_code, error_amount] = drawGaussianCurve3D( coeffs, sigma, imageIn, T, Idx, X, Y, Z)
        % Draw a gaussian curve in 3D using numerical integration

        if isempty( coeffs) 
            error('drawGaussianCurve3D: coeffs not provided')
        end
        if numel( size(imageIn) ) ~= 3 || length( coeffs) ~= 3 || length(sigma)~= 3
            error('drawGaussianCurve3D: all data must be 3-dimensional')
        end

        dim = length( size(imageIn) );
        error_code = 0;
        error_amount = 0;

        % Check if microtubule escapes the z-planes. Exit with error.
        %zends = polyval( coeffs{3}, T);
        %if any(zends <1) || any( zends > size(imageIn, 3))
            %imageCurve = imageIn;
            %error_code = 1;
            %return
        %end

        % Create coordinates for numerical integration
        % create parametric coord
        speedUp = 1.0;
        Tvec = T(1) : 0.01*speedUp : T(2);
        if length(Tvec) == 1
            Tvec = [Tvec(1) Tvec(1)+0.01];
        elseif length(Tvec) > 2
            Tvec( end) = [];
        end

        % Now for each value of the parameter, we'll do a gauss
        % quadrature numerical integration 
        DeltaT = median( diff(Tvec) );

        % get the offsets for gauss quadrature
        poff1 = ( (1/2) - sqrt(3)/6 ) * DeltaT;
        poff2 = ( (1/2) + sqrt(3)/6 ) * DeltaT;
        Tsort = sort([  T(1), Tvec + poff1,  Tvec + poff2, T(2)] );

        % Coefficients
        XCoef = coeffs{1}; YCoef = coeffs{2}; ZCoef = coeffs{3}; 
        xLoc = polyval( XCoef, Tsort )';
        yLoc = polyval( YCoef, Tsort )';
        zLoc = polyval( ZCoef, Tsort)'; 

        % re-parametrize by length
        dParam = 1;
        for jk = 1 : length( xLoc) -1
            dParam = [ dParam, dParam(end) + sqrt( diff( xLoc( jk:jk+1)).^2 + diff(yLoc(jk:jk+1)).^2 + diff(zLoc(jk:jk+1)).^2)]; 
        end
        dParamNew = linspace( dParam(1), dParam(end), length(dParam) );

        % fit once again, and find now parametric coordinates
        if any(isnan(dParam)) || any(isnan(xLoc))
            stoph = 1;
        end
        fitX = fit( dParam', xLoc, 'linearinterp');
        fitY = fit( dParam', yLoc, 'linearinterp');
        fitZ = fit( dParam', zLoc, 'linearinterp'); 
        xDat = fitX( dParamNew);
        yDat = fitY( dParamNew);
        zDat = fitZ( dParamNew); 

        if any(zDat < 1) || any(zDat >= size(imageIn,3))
            [xDat,yDat,zDat,error_amount] = TrimAndGetErrorZ( xDat,yDat,zDat,imageIn);
        end

        PtGaussLoc = [ xDat , yDat , zDat];
        imageCurve = NumericalConv( sigma, PtGaussLoc, {X,Y,Z}, imageIn, Idx );

        % TrimAndGetErrorZ {{{
        function [X,Y,Z,error_amount] = TrimAndGetErrorZ( X,Y,Z,imageIn)
            rmIdx = [];
            for jz = 1: length(Z)
                if Z(jz) < 1 || Z(jz) > size(imageIn,3)
                    rmIdx = [rmIdx, jz];
                end
            end
            error_amount = length(rmIdx)/length(X);
            X( rmIdx) = [];
            Y( rmIdx) = [];
            Z( rmIdx) = [];
        end
        % }}}

    end
    % }}}
    
    % drawGaussianCurve3D {{{
    function [imageCurve, error_code, error_amount] = drawGaussianCurve2DCoords( coords, sigma, imageIn, Idx, X, Y)
        % Draw a gaussian curve in 3D using numerical integration

        if isempty( coords)
            error('drawGaussianCurve2DCoords: coords not provided')
        end
        if numel( size(imageIn) ) ~= 2 || size( coords,1) ~= 2 || length(sigma)~= 2
            error('drawGaussianCurve2DCoords: all data must be 2-dimensional')
        end

        dim = length( size(imageIn) );
        error_code = 0;
        error_amount = 0;
        
        imageCurve = NumericalConv( sigma, [ coords(1,:)' , coords(2,:)'], {X,Y}, imageIn, Idx );

    end
    % }}}
    
    % drawGaussianCurve3D {{{
    function [imageCurve, error_code, error_amount] = drawGaussianCurve3DCoords( coords, sigma, imageIn, Idx, X, Y, Z)
        % Draw a gaussian curve in 3D using numerical integration

        if isempty( coords)
            error('drawGaussianCurve3DCoords: coords not provided')
        end
        if numel( size(imageIn) ) ~= 3 || size( coords,1) ~= 3 || length(sigma)~= 3
            error('drawGaussianCurve3D: all data must be 3-dimensional')
        end

        dim = length( size(imageIn) );
        error_code = 0;
        error_amount = 0;

        % Check if microtubule escapes the z-planes. Exit with error.
        %zends = polyval( coeffs{3}, T);
        %if any(zends <1) || any( zends > size(imageIn, 3))
            %imageCurve = imageIn;
            %error_code = 1;
            %return
        %end

        if any(coords(3,:) < 1) || any(coords(3,:) >= size(imageIn,3))
            [xDat,yDat,zDat,error_amount] = TrimAndGetErrorZ( coords(1,:), coords(2,:), coords(3,:),imageIn);
        else
            xDat = coords(1,:);
            yDat = coords(2,:);
            zDat = coords(3,:);
        end
        
        imageCurve = NumericalConv( sigma, [ xDat' , yDat' , zDat'], {X,Y,Z}, imageIn, Idx );

        % TrimAndGetErrorZ {{{
        function [X,Y,Z,error_amount] = TrimAndGetErrorZ( X,Y,Z,imageIn)
            rmIdx = [];
            for jz = 1: length(Z)
                if Z(jz) < 1 || Z(jz) > size(imageIn,3)
                    rmIdx = [rmIdx, jz];
                end
            end
            error_amount = length(rmIdx)/length(X);
            X( rmIdx) = [];
            Y( rmIdx) = [];
            Z( rmIdx) = [];
        end
        % }}}

    end
    % }}}

    % NumericalConv {{{
    function imConv = NumericalConv( Sig, pts, coordsCell, imPlane, idx )

%         p = gcp('nocreate');
%         if ~isempty( p)
%             nw = p.NumWorkers;
%             ip = round(linspace(1, length(idx), nw+1));
%             imConv = zeros([ size(imPlane), nw]);
%             
%             % Create copies of image volume
% %             imConv = imPlane .* 0;
%             fitdim = length(Sig); if length( coordsCell) ~= fitdim, error('Issue with passing grid indexes to NumericalConv2.m.'), end
%             NumGauss = size( pts, 1);
% 
%             CoordsX = coordsCell{1} .* ones( 1, NumGauss); CoordsY = coordsCell{2} .* ones( 1, NumGauss);
%             xx = pts( :, 1); yy = pts( :, 2);
%             s1 = Sig(1); s2 = Sig(2);
%             if fitdim==3,  
%                 CoordsZ = coordsCell{3} .* ones( 1, NumGauss);
%                 zz = pts( :, 3);s3 = Sig(3);
%             end
%             
%             parfor jw = 1: nw
%                 imT = 0*imPlane;
%                 % Create copies of image volume
%                 Xvals = exp( -( xx' - CoordsX(ip(jw):ip(jw+1),:) ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
%                 Yvals = exp( -( yy' - CoordsY(ip(jw):ip(jw+1),:) ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);
%                 if fitdim == 2
%                     imT( idx(ip(jw):ip(jw+1))') = sum( Xvals .* Yvals , 2);
%                 elseif fitdim == 3
%                     Zvals = exp( -( zz' - CoordsZ(ip(jw):ip(jw+1),:) ).^2 / (2 * s3^2) ) / ( 2 * s3^2 * pi)^(1/2);
%                     imT( idx(ip(jw):ip(jw+1))') = sum( Xvals .* Yvals .* Zvals , 2);
%                 end
%                 imConv(:,:,:,jw) = imT;
%             end
%             imConv = sum( imConv, 4);
%             
%         else
%                 
%             % Create copies of image volume
%             imConv = imPlane .* 0;
%             fitdim = length(Sig); if length( coordsCell) ~= fitdim, error('Issue with passing grid indexes to NumericalConv2.m.'), end
%             NumGauss = size( pts, 1);
% 
%             CoordsX = coordsCell{1} .* ones( 1, NumGauss); CoordsY = coordsCell{2} .* ones( 1, NumGauss);
%             xx = pts( :, 1); yy = pts( :, 2);
%             s1 = Sig(1); s2 = Sig(2);
%             if fitdim==3,  
%                 CoordsZ = coordsCell{3} .* ones( 1, NumGauss);
%                 zz = pts( :, 3);s3 = Sig(3);
%             end
%             Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
%             Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);
%             if fitdim == 2
%                 imConv( idx') = sum( Xvals .* Yvals , 2);
%             elseif fitdim == 3
%                 Zvals = exp( -( zz' - CoordsZ ).^2 / (2 * s3^2) ) / ( 2 * s3^2 * pi)^(1/2);
%                 imConv( idx') = sum( Xvals .* Yvals .* Zvals , 2);
%             end
%             
%         end
        
        % Create copies of image volume
        imConv = imPlane .* 0;

        fitdim = length(Sig); if length( coordsCell) ~= fitdim, error('Issue with passing grid indexes to NumericalConv2.m.'), end

        NumGauss = size( pts, 1);

        CoordsX = coordsCell{1} .* ones( 1, NumGauss);
        CoordsY = coordsCell{2} .* ones( 1, NumGauss);
        xx = pts( :, 1);
        yy = pts( :, 2);
        s1 = Sig(1);
        s2 = Sig(2);
        if fitdim==3,  
            CoordsZ = coordsCell{3} .* ones( 1, NumGauss);
            zz = pts( :, 3);
            s3 = Sig(3);
        end
        
        Xvals = exp( -( xx' - CoordsX ).^2 / (2 * s1^2) ) / ( 2 * s1^2 * pi)^(1/2);
        Yvals = exp( -( yy' - CoordsY ).^2 / (2 * s2^2) ) / ( 2 * s2^2 * pi)^(1/2);
        if fitdim == 2
            imConv( idx') = sum( Xvals .* Yvals , 2);
        elseif fitdim == 3
            Zvals = exp( -( zz' - CoordsZ ).^2 / (2 * s3^2) ) / ( 2 * s3^2 * pi)^(1/2);
            imConv( idx') = sum( Xvals .* Yvals .* Zvals , 2);
        end

    end
    % }}}

end
