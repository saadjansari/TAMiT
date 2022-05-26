function [origin, thetaInit, cf, L, cc] = get_tan_curvature_coeffs( coords, model_curvature, ax)
    % Parametric Polynomial representation
    % x(t) = a0 + a1*t + a2*t^2 + a3*t^3 + ... 
    % y(t) = b0 + b1*t + b2*t^2 + b3*t^3 + ... 
    % z(t) = c0 + c1*t + ... 
    
    % Order FIXED
    orderPoly = 4;
    
    % Input Parameter Checks
    
    % Check is axes is actually an axes object
    if nargin == 3
        try
            isAxes = strcmp(get(ax, 'type'), 'axes');
        catch
            isAxes = false;
        end
        if ~isAxes
            error('Input Value Error: ax must be of type axes')
        end
    end
   
    
    % Functions
    % Arc length parameter as cum distance between points.
    arc_length = @(x,y,z) cumsum( [0, sqrt( diff( x).^2 + diff( y).^2 + diff( z).^2 )]);
    % Length
    leng = @(x,y,z) sum( sqrt( diff( x).^2 + diff( y).^2 + diff( z).^2 ) );
    
    % Estimated coordinates
    X = coords(1,:);
    Y = coords(2,:);
    Z = coords(3,:);
    
    % Initial arc length (may be uneven spacing)
    t_uneven = arc_length(X,Y,Z);
    
    % Resampled coordinates with even spacing
    t_even = 0:0.25:round( max(t_uneven) );
    t_even_twice = 0:0.25:1.5*round( max(t_uneven) );
%     Xnew_twice = interp1( t_uneven, X, t_even_twice, 'spline','extrap' );
%     Ynew_twice = interp1( t_uneven, Y, t_even_twice, 'spline','extrap'  );
%     Znew_twice = interp1( t_uneven, Z, t_even_twice, 'spline','extrap'  );
    
    Xnew = interp1( t_uneven, X, t_even, 'spline','extrap' );
    Ynew = interp1( t_uneven, Y, t_even, 'spline','extrap'  );
    Znew = interp1( t_uneven, Z, t_even, 'spline','extrap'  );
    
    cc = [Xnew; Ynew; Znew];

    L = leng(Xnew,Ynew,Znew);
    
    % Get polynomial Coefficients
    cf1 = CurvedMT.estimatePolyCoefficients( [Xnew;Ynew;Znew], [orderPoly orderPoly 1], t_even);
    a0 = cf1{1}(end);
    a1 = cf1{1}(end-1);
    a2 = cf1{1}(end-2);
    a3 = cf1{1}(end-3);
    a4 = cf1{1}(end-4);
    b0 = cf1{2}(end);
    b1 = cf1{2}(end-1);
    b2 = cf1{2}(end-2);
    b3 = cf1{2}(end-3);
    b4 = cf1{2}(end-4);
    c0 = cf1{3}(end);
    c1 = cf1{3}(end-1);
    
    % Define Differentials
    x_p = a1 + 2*a2*t_even + 3*a3*(t_even.^2) + 4*a4*(t_even.^3);
    x_pp = 2*a2 + 6*a3*t_even + 12*a4*(t_even.^2);
    y_p = b1 + 2*b2*t_even + 3*b3*(t_even.^2)+ 4*b4*(t_even.^3);
    y_pp = 2*b2 + 6*b3*t_even + 12*b4*(t_even.^2);
    
    K = (x_p.*y_pp - y_p.*x_pp) ./ (x_p.^2 + y_p.^2 ).^(1.5);
    K_new = K(end)*ones( size(t_even_twice)); K_new(1:length(K)) = K;
    
    % Fourier fit the curvature function
    [f,g] = fit(t_even',K',model_curvature);
    [f,g] = fit(t_even_twice',K_new',model_curvature);
    if g.rsquare < 0.9
        warning('%s curvature fit might be bad!\nR^2 value for curvature model "%s" < 0.95 (R^2 = %.3d)',...
            model_curvature, model_curvature,g.rsquare)
    end
    cf = coeffvalues(f);
    
            
	% Get coordinates from coeffs (useful for plotting and looking at
	% accuracy of fitting) - These should match up well with the initial
	% coordinates
    x1 = polyval( cf1{1}, t_even); y1 = polyval( cf1{2}, t_even);
    
    % Get origin
    origin = coords(:,1);
    
    % Get initial tangent vector and initial theta vector (assume no
    % initial Z angle)
    tanInit = [a1, b1, c1];
    thetaInit = [atan2( tanInit(2), tanInit(1) ), pi/2];
    
    % Display
    if nargin == 3
        axes( ax); colormap gray; hold on;
        plot( coords(1,:), coords(2,:), 'r-', 'LineWidth',2, 'DisplayName', 'Detected Points')
        plot( x1,y1, 'g-', 'LineWidth',2, 'DisplayName', 'Polynomial Points')
        legend()
    end
end