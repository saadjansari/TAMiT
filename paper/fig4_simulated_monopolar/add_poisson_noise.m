 function imNoisy = add_poisson_noise(img, snr)
    noisyy = poissrnd( 5, size(img))/5;
    imNoisy = img + (1/snr)*noisyy;
end