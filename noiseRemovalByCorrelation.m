function imClean = noiseRemovalByCorrelation( imRaw, imMask3D)
% imRaw is a 3D image. The correlation is performed between the n-slices in the 3rd dimension. imClean is a 2D image for now (later on this could be made into a 3D image) by using it as a weight multiplier
if nargin == 1
    imMask2D = ones( size(imRaw(:,:,1)) );
else
    imMask2D = imMask3D(:,:, ceil(end/2) );
end

% for each layer of pixels in z, we take the brightest point as the source of the illumination, and compute theoretical values for the signal in every other z-slice. we then compute the correlation between the theoretical and the observed values. This returns a correlation score, which is a weighing used for noise reduction

psfz = 1.5;
numX = size(imRaw, 2);
numY = size(imRaw, 1);
numZ = size(imRaw, 3);
imCorr = zeros( numY, numX);
imClean = 0*imRaw;
% get list of x,y coordinates of pixels within the mask
idxYes = find(imMask2D); [y, x] = ind2sub( size(imMask2D), idxYes);
for j = 1 : length(x) 
    xx = x(j); yy = y(j);
    zVals = squeeze( imRaw(yy, xx, :) );
    [ampG, meanZ] =  max( zVals);
    coordZ = 1 : numZ;
    gaussTheory = exp( (coordZ - meanZ).^2 / psfz.^2 );
    gaussTheory = ( mat2gray( gaussTheory) * ampG)';
    
    % find correlation between expected and observed values
    imCorr( yy, xx) = corr2( zVals, gaussTheory);
end;
    imCorr = mat2gray(imCorr);
for jZ = 1 : numZ
    imClean( :,:,jZ) = imCorr .* imRaw(:,:,jZ);
end
dispImg( imCorr);
dispImg( max( imRaw, [], 3).*imMask2D, max(imClean, [], 3), [1 2])




end
