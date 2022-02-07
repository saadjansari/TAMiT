function save_as_tiff_stack(img3D,path)
% save_as_tiff_stack: Saves img3D as a tiff stack
%
% Input Parameters
% -----------------
%
% 1) img3D (3D array): this is the array that you wish to save as a tiff
% stack.
% 2) path (string): this is the location path to the tiff file. It must
% contain the full path to the file including the file name and a .tif
% extension.
%    E.g. path = '/path/to/correct/folder/image.tif'

    nSlices = size(img3D,3);
    imwrite(img3D(:,:,1),path)
    for jSlice = 2: nSlices
        imwrite(img3D(:,:,jSlice),path, 'WriteMode','append');
    end
end

