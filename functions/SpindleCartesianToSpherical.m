function [ len, phi, theta] = SpindleCartesianToSpherical( Aster1, Aster2 )
% SpindleCartesianToSpherical : converts spindle endpoints from cartesian
% to spherical coordinates. Takes Aster1 as the origin...
%   Detailed explanation goes here
    
    % coordinates for aster 1
    x0 = Aster1( 1);
    y0 = Aster1( 2);
    z0 = Aster1( 3);
    
    % coordinates for aster 2
    x1 = Aster2( 1);
    y1 = Aster2( 2);
    z1 = Aster2( 3);

    % get length
    len = norm( [ x0, y0, z0 ]- [x1, y1, z1] );
    
    % get phi
    phi = atan2( y1-y0, x1-x0);
    
    % get theta
    theta = atan2( norm([x1-x0, y1-y0]), z1-z0);
    
end

