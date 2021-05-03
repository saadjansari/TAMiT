function amplitude = findAmplitudeAlongCurveCoords( imageIn, coord)
    % used to measure intensity in the image at the given coordinates

    dim = length (size( imageIn) );

    if dim ==2
        idx = sub2ind( size(imageIn), coord(2, :), coord(1, :) );
    elseif dim ==3
        idx = sub2ind( size(imageIn), coord(2, :), coord(1, :), coord(3, :) );
    end
    amplitude = imageIn( idx);

end