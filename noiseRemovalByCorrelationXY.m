function imClean = noiseRemovalByCorrelation( imRaw, imMask3D)
% imRaw is a 3D image. The correlation is performed between the n-slices in the 3rd dimension. imClean is a 2D image for now (later on this could be made into a 3D image) by using it as a weight multiplier
if nargin == 1
    imMask2D = ones( size(imRaw(:,:,1)) );
else
    imMask2D = imMask3D(:,:, ceil(end/2) );
end

% for each layer of pixels in z, we take the brightest point as the source of the illumination, and compute theoretical values for the signal in every other z-slice. we then compute the correlation between the theoretical and the observed values. This returns a correlation score, which is a weighing used for noise reduction

psf = 1.2;
numX = size(imRaw, 2);
numY = size(imRaw, 1);
numZ = size(imRaw, 3);
imCorr = zeros( numY, numX, numZ);
imClean = 0*imRaw;

doThis = 0;
if doThis
    % get list of x,y coordinates of pixels within the mask
    idxYes = find(imMask2D); [y, x] = ind2sub( size(imMask2D), idxYes);
    imGauss = mat2gray( fspecial( 'gaussian', 3, psf) );
    for jZ = 1 : numZ
    for j = 1 : length(x) 
        xx = x(j); yy = y(j);
        imOrg = imRaw(yy-1:yy+1, xx-1:xx+1, jZ);
        ampG =  max( imOrg(:) );

        % find location of max values and treat that as a source of gaussian intensity
        [ampG, indG] = max( imOrg(:) );
        [indY, indX] = ind2sub( size(imOrg), indG);
        
        % create gaussian matric whose center is at (indY, indX)
        yList = 1 : size(imOrg, 1);
        xList = 1 : size(imOrg, 2);
        imTheory = 0*imOrg;
        for j1 = xList; for j2 = yList
            imTheory( j2, j1) = exp( (j2-indY)^2 / psf^2 ) .*  exp( (j1-indX)^2 / psf^2 );
        end; end
        imTheory = mat2gray( imTheory) .* ampG;
        
        % find correlation between expected and observed values
        imCorr( yy, xx, jZ) = corr2( imOrg, imTheory);
    end;
    end;

    imCorr = mean( imCorr, 3);
    imCorr = imCorr + abs( min( imCorr(:) ) ); 
    imCorr = mat2gray(imCorr .* imMask2D);
    for jZ = 1 : numZ
        imClean( :,:,jZ) = imCorr .* imRaw(:,:,jZ);
    end
    dispImg( max( imRaw, [], 3).*imMask2D, max(imClean, [], 3), [1 2])
end


mu =  mean(imRaw(imMask3D) );
sig = std(imRaw(imMask3D) ); 

PSF = fspecial('gaussian', 8, 1.3);
imPlane = mean( imRaw, 3);
imPlane( imPlane == 0) = median( imPlane( imPlane ~= 0) );
imClean = deconvlucy( imPlane, PSF, 45, 3*sig).*imMask2D;
% dispImg( imPlane.*imMask2D, imClean, [1 2])
% figure; surf( [imPlane.*imMask2D , imClean]); set(gcf, 'WindowState', 'maximized');

end
