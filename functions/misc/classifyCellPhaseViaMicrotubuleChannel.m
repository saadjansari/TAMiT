function cellphase = classifyCellPhaseViaMicrotubuleChannel( mtChannel, times)
	% cell phase classification will be based on the tubulin channel

    size(mtChannel)

    cellphase = cell( size(times));


	% we care about the distribution of intensity across pixels inside the cell
    figure;
    for jT = times

        imXYZ = mtChannel(:,:,:, jT);
        imXY = mean( imXYZ, 3);
	    imXY( imXY == 0) = median( imXY( imXY~= 0) );
        imXY = imgaussfilt( imXY, 1);

        % find SNR
        pixVal = imXY(:);
        noise = median( pixVal );
        snr = pixVal/noise; 
        numHighSNR = sum(snr > 3);
        max_snr = max( snr);
        
        if numHighSNR > 10
            cellphase{jT} = 'Met';
        else cellphase{jT} = 'Int'; 
        end

%         histogram( imXY(:) ); set(gca, 'yscale', 'log')

        imagesc( imXY); axis equal; colormap gray
        title( sprintf( 't = %d , maxSNR = %.2f , numHighSNR = %d', jT, max_snr, numHighSNR) )
        pause(0.2)

    end

end 
