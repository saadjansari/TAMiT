function [origin, thetaInit, nV, L, cX1, cY1, cZ1] = get_tan_normal_curve_coords( coords, orderXY, ax)
    % Parametric Polynomial representation
    % x(t) = a0 + a1*t + a2*t^2 + a3*t^3 + ... 
    % y(t) = b0 + b1*t + b2*t^2 + b3*t^3 + ... 
    % z(t) = c0 + c1*t + ... 
    
    % Input Parameter Checks
    % Check value of order
    if orderXY ~= 2 && orderXY ~= 3
        error('Input Value Error: orderXY is not 2 or 3')
    end
    
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
    t_even = 0:round( max(t_uneven) );
    Xnew = interp1( t_uneven, X, t_even );
    Ynew = interp1( t_uneven, Y, t_even );
    Znew = interp1( t_uneven, Z, t_even );
%     [Xnew, t_even] = resample(X,t_uneven,1);
%     [Ynew, t_even] = resample(Y,t_uneven,1);
%     [Znew, t_even] = resample(Z,t_uneven,1);
    L = leng(Xnew,Ynew,Znew);
    cX1 = Xnew; cY1 = Ynew; cZ1 = Znew;
    
    % Get length
%     L = sum( sqrt( diff( coords(1,:)).^2 + diff( coords(2,:)).^2 + diff( coords(3,:)).^2 ) );
    
    % Average spacing between coords
%     meanSpacing = ceil( L/ length( coords(1,:) ) );
    
    % Interpolate coordinates finely (meanSpacing number of points are
    % added between consecutive coordinates)
%     [cX1,cY1,cZ1] = Methods.InterpolateCoords3( coords(1,:), coords(2,:), coords(3,:), meanSpacing );
    
    
    % Get polynomial Coefficients
%     cf1 = CurvedMT.estimatePolyCoefficients( [cX1;cY1;cZ1], [orderXY orderXY 1], t2);
    cf1 = CurvedMT.estimatePolyCoefficients( [Xnew;Ynew;Znew], [4 4 1], t_even);
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
%     a3 = 0;
%     b3 = 0;
    
    % Define Differentials
    x_p = a1 + 2*a2*t_even + 3*a3*(t_even.^2) + 4*a4*(t_even.^3);
    x_pp = 2*a2 + 6*a3*t_even + 12*a4*(t_even.^2);
    y_p = b1 + 2*b2*t_even + 3*b3*(t_even.^2)+ 4*b4*(t_even.^3);
    y_pp = 2*b2 + 6*b3*t_even + 12*b4*(t_even.^2);
    
    K = (x_p.*y_pp - y_p.*x_pp) ./ (x_p.^2 + y_p.^2 ).^(1.5);
    
    % Fourier fit the curvature function
    model_curvature = 'fourier2';
    [f,g] = fit(t_even',K',model_curvature);
    if g.rsquare < 0.95
        warning('%s curvature fit might be bad!\nR^2 value for curvature model "%s" < 0.95 (R^2 = %.3d)',...
            model_curvature, model_curvature,g.rsquare)
    end
    
    
            
	% Get coordinates from coeffs (useful for plotting and looking at
	% accuracy of fitting) - These should match up well with the initial
	% coordinates
%     x1 = polyval( cf1{1}, t2); y1 = polyval( cf1{2}, t2);
    x1 = polyval( cf1{1}, t_even); y1 = polyval( cf1{2}, t_even);
    
    % Get origin
    origin = coords(:,1);
    
    % Get initial tangent vector and initial theta vector (assume no
    % initial Z angle)
    tanInit = [a1, b1, c1];
    %thetaInit = [atan2( tanInit(2), tanInit(1) ), pi/2];
    
    % Normal vector coefficients
%     nV = 1;
%     switch orderXY
%         case 2
%             
%             % There is a single coefficient
%             nV = 2*(a1*b2 - a2*b1);
%             
%         case 3
%             a3 = cf1{1}(end-3);
%             b3 = cf1{2}(end-3);
%             
%             nV = [ 2*(a1*b2 - a2*b1), ...
%                     6*(a1*b3 - a3*b1),...
%                     6*(a2*b3 - a3*b2)];
%             
%         otherwise
% %             error('orderXY was neither 2 nor 3')
%     end
    
    % Display
    if nargin == 3
        axes( ax); colormap gray; hold on;
        plot( coords(1,:), coords(2,:), 'r-', 'LineWidth',2, 'DisplayName', 'Detected Points')
        plot( x1,y1, 'g-', 'LineWidth',2, 'DisplayName', 'Polynomial Points')
        legend()
    end
end