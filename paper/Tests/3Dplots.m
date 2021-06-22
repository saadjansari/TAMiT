
% XZ plot
imXZ = squeeze(max(Image2Fit, [],2));
intMin = min( Image2Fit(:) ); intMax = max( Image2Fit(:) );
nX = size(imXZ,2);
nY = size(imXZ,1);

imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

figure;
ax = tight_subplot(1, 2, -0.7);
drawnow; pause(1)
axes( ax(1) );
imagesc( imXZ ); eval( imageSets);
axes( ax(2) );
imagesc( imXZ ); eval( imageSets); hold on;
fitEngine.GetFeature().displayFeatureXZ(gca);
title('XZ Features')

% YZ plot
imYZ = squeeze(max(Image2Fit, [],1));
intMin = min( Image2Fit(:) ); intMax = max( Image2Fit(:) );
nX = size(imYZ,2);
nY = size(imYZ,1);

imageSets = 'colormap gray; axis equal; axis ij; set( gca, ''xlim'', [1 nX], ''ylim'', [1 nY], ''XTick'', [], ''YTick'', [], ''CLim'', [intMin intMax], ''FontSize'', 14)';

figure;
ax = tight_subplot(1, 2, -0.7);
drawnow; pause(1)
axes( ax(1) );
imagesc( imaZ ); eval( imageSets);
axes( ax(2) );
imagesc( imYZ ); eval( imageSets); hold on;
fitEngine.GetFeature().displayFeatureYZ(gca);
title('YZ Features')
