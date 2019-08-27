% This is a script that plays a movie from a pre-segmented cell
%
% Plays the image movie of the Cell.
% cChannels can either be a scalar or a vector of upto 3 elements specifying the channels to play
% cColor can either be a string or a string array of upto 3 elements specifying the color of each of the channels. the size of cColors must match the size of cChannels

% load the movie from the segmented cell
cellpath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/1095_50msG_100msR_7Z_007_cells/1095_50msG_100msR_7Z_007_26.mat';
cellData = importSingleCell( cellpath );
image = cellData.cell3D;

numChannels = size(image, 5);
cChannels = 2;
% cChannels = 1 : numChannels;
colors = {'G', 'R', 'B'};
cColors = colors( cChannels);

if length( cChannels) > 3
    error('playMovie: cChannels cannot contain more than 3 elements')
end

img = max( image, [], 3);
img = squeeze( img);
img = permute( img, [ 1 2 4 3] ); % XYCT

imStack = zeros( size(img, 1), size(img,2), 3, size(img, 4) );
rgb = {'R', 'G', 'B'};

% Find cChannel for red color
cR = cChannels( find( strcmp(cColors, 'R') ) );
if ~isempty( cR)
    imStack(:,:,1,:) = mat2gray( img(:,:,cR, :) );
end

% Find cChannel for green color
cG = cChannels( find( strcmp(cColors, 'G') ) );
if ~isempty( cG)
    imStack(:,:,2,:) = mat2gray( img(:,:,cG, :) );
end

% Find cChannel for blue color
cB = cChannels( find( strcmp(cColors, 'B') ) );
if ~isempty( cB)
    imStack(:,:,3,:) = mat2gray( img(:,:,cB, :) );
end

fps = 10;
implay( imStack, fps);

