function save_as_tiff_stack(img3D,path)
%save_as_tiff_stack: Saves img3D as a tiff stack

    nSlices = size(img3D,3);
    imwrite(img3D(:,:,1),path)
    for jSlice = 2: nSlices
        imwrite(img3D(:,:,jSlice),path, 'WriteMode','append');
    end
end

