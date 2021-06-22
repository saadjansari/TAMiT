function masks = find_nuclear_mask( frames)

    masks = cell(1, length(frames));
    for jf = 1 : length(frames)
        Image2Fit = frames{jf};
        im3 = imgaussfilt3( Image2Fit, 1);
        im2 = mean( im3,3);
        
        im2 = mat2gray(im2);

        func = 'histeq';
        % HISTEQ
        imEdit = feval( func, im2);

        %Threshold
        thresh = multithresh(imEdit,2);
        imBinary = imbinarize( imEdit, thresh(1) );

        % SE
        se1 = strel('disk', 1);
        se3 = strel('disk', 3);
        im_dilated = imdilate( imerode( imBinary, se3), se3);
        im_filled = imfill( im_dilated, 'holes');
        im_filled = mat2gray(im_filled);
        imColor = cat( 3, im_filled, zeros( size(im2) ), im2);

        %montage({im2,imEdit,imBinary,im_dilated, im_filled, imColor},'Size',[1 6])
        % title("Original Image and Enhanced Images using imadjust, histeq, and adapthisteq")
        mask3D = zeros( size(Image2Fit) );
        for jz = 1 : size(mask3D,3)
            mask3D(:,:,jz) = im_filled;
        end
        masks{jf} = mask3D;
    end
end
