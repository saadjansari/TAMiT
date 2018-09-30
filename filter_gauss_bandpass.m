function imfilter = filter_gauss_bandpass( image, sigMax, sigMin)
% uses a bandpass filter to filter the image

    if iscell(image)
        cellYes = 1;
        imfilter = cell( size(image) );
        image = cat( 3, image{:});
    else, 
        cellYes = 0; 
        imfilter = 0*image;
    end

    % for each slice, do a 2d bandpass
    for jT = 1 : size(image,3)

        im_filter_highfreq = imgaussfilt( image(:,:,jT), sigMin);
        im_filter_lowfreq = imgaussfilt( image(:,:,jT), sigMax);
        imgFilt = mat2gray(im_filter_highfreq - im_filter_lowfreq);
        if cellYes, imfilter{jT} = imgFilt;  else, imfilter(:,:,jT) = imgFilt; end
%         dispImg( image(:,:,jT), imgFilt, [1 2]);

    end

end
