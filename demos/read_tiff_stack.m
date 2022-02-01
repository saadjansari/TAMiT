function tiff_stack = read_tiff_stack(path)
%read_tiff_stack: Read an image stack from tiff

    tiff_info = imfinfo(path); % return tiff structure, one element per image
    tiff_stack = imread(path, 1) ; % read in first image
    
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread(path, ii);
        tiff_stack = cat(3 , tiff_stack, temp_tiff);
    end
end