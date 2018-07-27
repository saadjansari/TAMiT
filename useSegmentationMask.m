function IsolatedCells = useSegmentationMask( image2D, imageMask, method)
% useSegmentationMask : applies a segmentation mask on an image and returns
% an image of a specific cell/cells;
%
%   method - 'prompt' prompts the user to select cells
%          - [ 1, 5, 7, ...] selects the cells within this vector
%
%   imageMask can either be a logical Mask or a labeled Mask
%
%
%
%
%   Detailed explanation goes here

% create logical and labeled segmentation masks. Find number of cells in
% this mask
imLogic = logical( imageMask);
imLabel = bwlabel( imLogic);
numCells = max( imLabel(:) );

% if image2D is a 3D image with the third dimension being time, then add up
% the times
numTimes = size(image2D, 3);
if numTimes > 1
    image3D = image2D;
    image2D = mat2gray( sum( image2D, 3) );
elseif numTimes == 1
    image3D = image2D;
end

% We can go one of two routes. If the method is 'prompt', then we'll
% display an image for the user and ask the user to select cells to
% isolate. Otherwise, we will just isolate the cells for the user that
% specified in the method vector.

 % Create the image with numbers on top of cells to prompt the user
% with. To do this, we need the centroid of each cell which which
% dictate the position the number will be placed on the image.
cc = bwconncomp( imLogic);
stats = regionprops( cc, 'Centroid');

% create background image for user
imUser = image2D .* imLogic;
idx = (1:numCells)';
pos = zeros( numCells, 2);
for jObj = 1 : numCells
    pos(jObj, :) = round( stats(jObj).Centroid );
end

imUser( imUser > 0.3) = 0.3; imUser = mat2gray(imUser);
T = multithresh( imUser( imUser>0), 2);
%     imUser( imUser < T(1)) = 0; imUser = mat2gray(imUser);
% add numbers to the image
imUserNumbered = insertText( imUser, pos, idx, 'AnchorPoint', 'Center', 'BoxOpacity', 0, 'TextColor', 'red', 'FontSize', 20);


if strcmp( method, 'prompt')
       
    % display the image and prompt the user to select cells
    h = figure;
    imagesc( imUserNumbered)
    title('Cell Selection Window')
    set( h, 'WindowStyle', 'Docked');

    % Prompt user to select cells, and check that specified cells are
    % possible.
    flagInput = 1; count = 0;
    while flagInput
        CellsSelected = input(['\n Please specify a single cell to select for analysis : \n'...
            '    Valid cell numbers are ' num2str( 1) '-' num2str( numCells ) ...
            '. \n    Selected Cells = ']);
        count = count+1;
        if all( ismember( CellsSelected, [1 : numCells] ) )
           flagInput = 0;
           close(h)
        elseif count > 5
            close(h)
            error('useSegmentationMask: Too many user tries, software aborting. \n Please run it again. \n')
        else
           fprintf('\n Incorrect Entry! Please try again. Please specify cells in the provided range and in this format [1, 5, 24, ...] \n') 
        end
    end

elseif isa( method, 'double') && size(method, 1)==1
    
    CellsSelected = method;
    
else
    error('useSegmentationMask: Incorrect argument, method must be either "prompt" or a vector specifying the cells to isolate')
end

% Sort the Cells selected so that they are in order of number
CellsSelected = sort(CellsSelected, 'ascend');

% Now we have selected some cells. For each cell we'll crop a 150*150 box
% centered at the cell.
boxSize = 150;
numWanted = length( CellsSelected);

% We need the centroid of the cells to crop a region right around it.
% We will also ensure that the isolated cells are elliptical
imLogicPad = padarray( imLogic, [boxSize, boxSize], 0);
cc = bwconncomp( imLogicPad);
stats = regionprops( cc, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');

% pad the original image
image3Dpad = padarray( image3D, [boxSize, boxSize],0 );
    
for jCell = 1 : numWanted
    
    % Turn on only the correct cell label
    currCell = CellsSelected( jCell);
    imMask = imLabel;
    imMask( imMask ~= currCell ) = 0;
    imMask = logical( imMask);
    
    % Create elliptical mask
    t = linspace(0,2*pi,1000);
    a = stats( currCell).MajorAxisLength/2;
    b = stats( currCell).MinorAxisLength/2;
    Xc = stats( currCell).Centroid(1);
    Yc = stats( currCell).Centroid(2);
    phi = deg2rad( -stats( currCell).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    idx = sub2ind( size( imLogicPad), round(y), round(x) );
    imMaskPad = 0*imLogicPad;
    imMaskPad( idx) = 1; imMaskPad = imfill( imMaskPad, 'holes');
    
    % dilate it just a little to ensure any edge features are captured.
    imMaskPad = imdilate( imMaskPad, strel('disk', 1) );
    
    % multiply with the original image to recover intensity information
    imMask3DPad = repmat( imMaskPad, 1,1, numTimes);
    imCell3DPad = imMask3DPad .* mat2gray( image3Dpad);
    
    % get region around cell;
    ub = round( stats( currCell).Centroid + boxSize/2 );
    lb = round( stats( currCell).Centroid - boxSize/2 );
    
    % obtain cropped regions for both the segmentation version and the
    % field of view region
    IsolatedCells( jCell).cell3D = imCell3DPad( lb(2): ub(2), lb(1): ub(1), : );
    IsolatedCells( jCell).cellMIP = max( IsolatedCells( jCell).cell3D, [], 3);
    IsolatedCells( jCell).raw = mat2gray( image3Dpad( lb(2): ub(2), lb(1): ub(1), : ) );
    IsolatedCells( jCell).rawMIP = mat2gray( max( IsolatedCells( jCell).raw, [], 3) );
    IsolatedCells( jCell).cellNumber = currCell;
    IsolatedCells( jCell).locations = imUserNumbered;
    
end



end

