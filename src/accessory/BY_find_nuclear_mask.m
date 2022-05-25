function mask = BY_find_nuclear_mask( img3)

    % Check if image is 3D or 2D
    if length( size(img3) ) == 3
        dim = 3;
    else
        dim=2;
    end
    
    if dim==3
        img = imgaussfilt3( img3, 1);
    else
        img = imgaussfilt( img3,1);
    end
    img = mat2gray(img);
    
    % Define function for contrast adjustsment
    func = 'histeq'; % imadjust,histeq or adapthisteq
    
    % HISTEQ
    imEdit = feval( func, img);

    %Threshold
    thresh = multithresh(imEdit(:),2);
    imBinary = imbinarize( imEdit, thresh(1) );

    % SE
    se3 = strel('disk', 3);
    im_filled = zeros( size(imBinary) );
    for jZ = 1 : size(img3,3)
        im_dilated = imdilate( imerode( imBinary(:,:,jZ), se3), se3);
        im_filled(:,:,jZ) = imfill( im_dilated, 'holes');
    end
    im_filled = mat2gray(im_filled);

    %montage({im2,imEdit,imBinary,im_dilated, im_filled, imColor},'Size',[1 6])
    
    mask = im_filled;
end
