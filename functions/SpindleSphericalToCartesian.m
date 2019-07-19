function [ X, Y, Z] = SpindleSphericalToCartesian( len, phi, theta )
% SpindleSphericalToCartesian : converts spindle endpoints from spherical
% to cartesian coordinates. Takes Aster1 to be at the origin...
%   Detailed explanation goes here
    
    % coordinates for aster 1 (at origin)
    x0 = 0;
    y0 = 0;
    z0 = 0;

    % get coordinates for aster 2
    X = len .* cos( phi) .* sin(theta);
    Y = len .* sin( phi) .* sin(theta);
    Z = len .* cos( theta);

end

